: Ih current
: Created 8/6/02 - nwg
: for the formulation of Angelo's data, the original parameter correspond to a maximum value of 250 ms. In their supplement data, this value is about 325 ms. So I corrected a value1.3
NEURON {
	SUFFIX hpkj
	NONSPECIFIC_CURRENT i
	RANGE ghbar, eh, i
	GLOBAL ninf, ntau
:    THREADSAFE	
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}
CONSTANT {
	q10 = 3

}
PARAMETER {
	v	 	(mV)
	celsius (degC)
	ghbar = .0001	(S/cm2)

	eh = -30	(mV)
}

ASSIGNED {
	i (mA/cm2)
	qt
	ninf
	ntau
}

STATE {
	n
}

INITIAL {
    qt = q10^((celsius-22 (degC))/10 (degC))
	rates(v)
	n = ninf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	i = ghbar*n*(v - eh)
}

DERIVATIVE states {
	rates(v)
	n' = (ninf - n)/ntau
}

PROCEDURE rates(v (mV)) {
:	ninf = 1/(1+exp((v+90.3+3)/9.9))
	ninf = 1/(1+exp((v+90.3+3+3)/9.67))

	ntau = 1000/(0.62*(exp((v+68)/-22)+exp((v+68)/7.14)))/qt/1.3
}