: Ca-dependent K channels (BK) - alphabeta4 and alpha
: Bin Wang, Robert Brenner, and David Jaffe - Originally written May 27, 2010
: 
: June 1, 2010 - added double exponential function for voltage-dependent activation 
:
: July 3, 2010 - changed voltage-dependence for the two channels based on revised data
:
: April 2, 2011 - adjusted parameters based on updated Bin data
: modified by Yunliang October 19th, 2015
: Notice that the experiments are done under room temperature after checking with David Jaffe!!!!
: However, to make slow component function during the ISI, here the Q10 correction is not included.
NEURON {
	SUFFIX abBK
	USEION ca READ cai
	USEION k READ ek WRITE ik
	RANGE gabkbar,gabk, ik
:	THREADSAFE
}

UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)

}
CONSTANT {
    q10 = 3
}

PARAMETER {
	gabkbar = .01	(S/cm2)	: maximum permeability - alphabeta
    cai (mM)
    base = 1  	(mV)	: alphabeta4 base time constant
}

ASSIGNED {
	v		(mV)
	ek		(mV)
	ik		(mA/cm2)
    gabk		(S/cm2)
    abinf		(mV)
    abtau		(ms)
    qt
}

STATE { ab }

BREAKPOINT {
	SOLVE state METHOD cnexp
	gabk = gabkbar*ab
	ik = (gabk)*(v - ek)
}

DERIVATIVE state {	: exact when v held constant; integrates over dt step
	rates(v, cai)				      
	ab' = (abinf-ab)/abtau
}

INITIAL {
    qt = q10^((celsius-34 (degC))/10 (degC))
	rates(v, cai)
	ab = abinf
}

: alpha-beta4 channel properties


FUNCTION shiftab(cai (mM))  {
	shiftab = 25 - 55.7 + 136.9*exp(-.28*cai*1e3)
}


FUNCTION peakab(cai (mM))  {
	peakab = 13.7 + 234*exp(-.72*cai*1e3)
}

: Double sigmoid function for tau voltage-dependence


FUNCTION taufunc(v (mV)) {
	 taufunc = 1 / (          (10*(exp(-v/63.6) + exp (-(150-v)/63.6)))  - 5.2                  )
	 if (taufunc <= 0.2) {	  : stop the function between 0.2 and 1
	    taufunc = 0.2
	 }

}

PROCEDURE rates(v (mV), cai (mM)) {
	  LOCAL range, vv

	  : alpha-beta4 model

	  abinf = -56.449 + 104.52*exp(-.22964*cai*1e3) + 295.68*exp(-2.1571*cai*1e3)

	  abinf = 1/(1+exp((abinf-v)/(25/1.6)))

	  vv = v + 100 - shiftab(cai)
	  abtau = taufunc(vv)
	  range = peakab(cai)-base
	  abtau = ((range*((abtau-.2)/.8)) + base)/qt*2*2

}
