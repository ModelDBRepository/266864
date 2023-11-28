TITLE Voltage-gated low threshold potassium current from Kv1 subunits

COMMENT
The time kinetics are from the Original paper (Akeman reduce the time constatant to about 1/3);
activation curve are from original paper; inactivation data are from Stephane's paper Fig 6E 
(This current should be a combination of kv1.2 which is not sensitive to TEA and kv4.3 which is not completely blocked).
Not sure this steady state inactivation curve is right, but it would not affect the normal firing of PC. It is used to 
reproduce history dependence of development of secondary spikes. Kv4 seems not to be able to do it according to its fast inactivation kinetics.
Ik1 inactivates very slow. We set the inactivation time constant to 1000 ms. 
NEURON implementation of a potassium channel from Kv1.1 subunits
We should also refer to Gating, modulation and subunit composition of voltage-gated K+channels in dendritic inhibitory interneurones of rathippocampus
In this paper,they observe a slow delayed rectifier K current. This current is surprisingly not sensitive to 4AP,but sensitive to high does
TEA. Compared with the paper by Marco Martin 5698 • The Journal of Neuroscience, July 2, 2003 • 23(13):5698 –5707, the K currents measured in the
dendrite, about 20% of the K currents is still resistant to even 3 mM 4-AP. However, 10 mM TEA can decrease the K currents to only 10%, suggesting a K
current sensitive  to high concentration TEA instead of 4AP, similar with the work by CC Lien. I guess in the soma, there should be IK1.

Kinetic data taken from: Zerr et al., J.Neurosci. 18 (1998) 2842
Vhalf = -28.8 +/- 2.3 mV; k = 8.1 +/- 0.9 mV

The voltage dependency of the rate constants was approximated by:

alpha = ca * exp(-(v+cva)/cka)
beta = cb * exp(-(v+cvb)/ckb)

Parameters ca, cva, cka, cb, cvb, ckb
are defined in the CONSTANT block.

Laboratory for Neuronal Circuit Dynamics
RIKEN Brain Science Institute, Wako City, Japan
http://www.neurodynamics.brain.riken.jp

Reference: Akemann and Knoepfel, J.Neurosci. 26 (2006) 4602
Date of Implementation: April 2005
Contact: akemann@brain.riken.jp

ENDCOMMENT

NEURON {
	SUFFIX Kv1
	USEION k READ ek WRITE ik
	RANGE gk, gbar, ik
	GLOBAL ninf, taun
:    THREADSAFE
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(nA) = (nanoamp)
	(pA) = (picoamp)
	(S)  = (siemens)
	(nS) = (nanosiemens)
	(pS) = (picosiemens)
	(um) = (micron)
	(molar) = (1/liter)
	(mM) = (millimolar)		
}

CONSTANT {
	q10 = 3

	ca = 0.12889 (1/ms)
	cva = 45 (mV)
	cka = -33.90877 (mV)

	cb = 0.12889 (1/ms)
      cvb = 45 (mV)
	ckb = 12.42101 (mV)         
}

PARAMETER {
	v (mV)
	celsius (degC)
	
	gbar = 0.011 (mho/cm2)   <0,1e9>
}


ASSIGNED {
 	ik (mA/cm2) 
	ek (mV)
	gk  (mho/cm2)
	ninf
	hinf
	taun (ms)
	tauh
	alphan (1/ms)
	betan (1/ms)
	qt
}

STATE { n h}

INITIAL {
	qt = q10^((celsius-22 (degC))/10 (degC))
	rates(v)
	n = ninf
	h = hinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
      gk = gbar * n^4 *h
	ik = gk * (v - ek)
}

DERIVATIVE states {
	rates(v)
	n' = (ninf-n)/taun
	h' = (hinf-h)/tauh 
}

PROCEDURE rates(v (mV)) {
	alphan = alphanfkt(v)
	betan = betanfkt(v)
	ninf = alphan/(alphan+betan) 
	taun = 1/(qt*(alphan + betan))
:	hinf = 0.1765+0.8235/(1+exp((v+70)/11.5))
	tauh = 1000/qt
	hinf = 1/(1+exp((v+66.16)/6.1881))
:	tauh = 1000/(1+exp((v+58.72)/3.005))/qt
}

FUNCTION alphanfkt(v (mV)) (1/ms) {
	alphanfkt = ca * exp(-(v+cva)/cka) 
}

FUNCTION betanfkt(v (mV)) (1/ms) {
	betanfkt = cb * exp(-(v+cvb)/ckb)
}




