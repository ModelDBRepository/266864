TITLE Voltage-gated potassium channel from Kv3 subunits
COMMENT
Voltage-gated potassium channel with high threshold and fast activation/deactivation kinetics

KINETIC SCHEME: Hodgkin-Huxley (n^4)
n'= alpha * (1-n) - betha * n
g(v) = gbar * n^4 * ( v-ek )

The rate constants of activation (alpha) and deactivation (beta) were approximated by:

alpha(v) = ca * exp(-(v+cva)/cka)
beta(v) = cb * exp(-(v+cvb)/ckb)

Parameters can, cvan, ckan, cbn, cvbn, ckbn are given in the CONSTANT block.
Values derive from least-square fits to experimental data of G/Gmax(v) and taun(v) in Martina et al. J Neurophysiol 97:563-571, 2007
Model includes a calculation of Kv gating current

Reference: Akemann et al., Biophys. J. (2009) 96: 3959-3976
Notice that there is another set of data related with Kv3 by McKay and Turner European Journal of Neuroscience, Vol. 20, pp. 729–739, 2004
in that paper, the activation threshold of Kv3 is much lower.
Laboratory for Neuronal Circuit Dynamics
RIKEN Brain Science Institute, Wako City, Japan
http://www.neurodynamics.brain.riken.jp

Date of Implementation: April 2007
Contact: akemann@brain.riken.jp

ENDCOMMENT


NEURON {
	SUFFIX Kv3
	USEION k READ ek WRITE ik
	RANGE gbar, g, ik,vshift
	GLOBAL ninf, tau
:    THREADSAFE
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(nA) = (nanoamp)
	(pA) = (picoamp)
	(S)  = (siemens)
	(mS) = (millisiemens)
	(nS) = (nanosiemens)
	(pS) = (picosiemens)
	(um) = (micron)
	(molar) = (1/liter)
	(mM) = (millimolar)		
}

CONSTANT {
	e0 = 1.60217646e-19 (coulombs)
	q10 = 2.7

	ca = 0.22 (1/ms)
	cva = 16 (mV)
	cka = -26.5 (mV)
	cb = 0.22 (1/ms)
	cvb = 16 (mV)
	ckb = 26.5 (mV)
	
	zn = 1.9196 (1)		: valence of n-gate
}

PARAMETER {
	vshift = 0
	gbar = 0.005 (S/cm2)   <0,1e9>
}

ASSIGNED {
	celsius (degC)
	v (mV)
	
	ik (mA/cm2)
 
	ek (mV)
	g (S/cm2)
	qt (1)

	ninf (1)
	tau (ms)
	alpha (1/ms)
	beta (1/ms)
}

STATE { n }

INITIAL {
	qt = q10^((celsius-22 (degC))/10 (degC))
	rateConst(v)
	n = ninf
}

BREAKPOINT {
	SOLVE state METHOD cnexp
      g = gbar * n^4 
	ik = g * (v - ek)

}

DERIVATIVE state {
	rateConst(v)
	n' = alpha * (1-n) - beta * n
}

PROCEDURE rateConst(v (mV)) {
	alpha = qt * alphaFkt(v)
	beta = qt * betaFkt(v)
	ninf = alpha / (alpha + beta) 
	tau = 1 / (alpha + beta)
}

FUNCTION alphaFkt(v (mV)) (1/ms) {
	alphaFkt = ca * exp(-(v+cva+vshift)/cka)
}

FUNCTION betaFkt(v (mV)) (1/ms) {
	betaFkt = cb * exp(-(v+cvb+vshift)/ckb)
}






