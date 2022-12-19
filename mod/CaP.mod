TITLE P-type calcium channel

COMMENT

According to Benton&Raman data
lower threshold but relatively large time constant compared with Sungho's model (According to Bruce Bean)
Also the ssa is steep. In this model, it is better not to shift the SSA to the left.
time speeded up by 2 times May 9 2016 (no longer)
ENDCOMMENT

NEURON {
	SUFFIX newCaP
	USEION ca READ cai, cao WRITE ica
	RANGE pcabar, ica,vshift,kt
	GLOBAL minf, taum
	GLOBAL monovalConc, monovalPerm
:	THREADSAFE
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
	F = 9.6485e4 (coulombs)
	R = 8.3145 (joule/kelvin)

:	cv = 19 (mV)
:	ck = 5.5 (mV)
    cv = 30.5 (mV)
    ck = 4.113 (mV)
}

PARAMETER {
	v (mV)
	celsius (degC)

	cai (mM)
	cao (mM)
    vshift =0
	pcabar = 6e-5 (cm/s)
	monovalConc = 140 (mM)
	monovalPerm = 0
	kt=1
}

ASSIGNED {
	qt
	ica (mA/cm2)
      minf 
	taum (ms)
	T (kelvin)
	E (volt)
	zeta
}

STATE { m }

INITIAL {
	qt = q10^((celsius-22 (degC))/10 (degC))
	T = kelvinfkt( celsius )
	rates(v)
	m = minf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ica = (1e3) * pcabar * m * ghk(v, cai, cao, 2)
}

DERIVATIVE states {
	rates(v)
	m' = (minf-m)/taum
}

FUNCTION ghk( v (mV), ci (mM), co (mM), z )  (coulombs/cm3) { 
	E = (1e-3) * v
      zeta = (z*F*E)/(R*T)	
	
	: ci = ci + (monovalPerm) * (monovalConc) :Monovalent permeability

	if ( fabs(1-exp(-zeta)) < 1e-6 ) {
	ghk = (1e-6) * (z*F) * (ci - co*exp(-zeta)) * (1 + zeta/2)
	} else {
	ghk = (1e-6) * (z*zeta*F) * (ci - co*exp(-zeta)) / (1-exp(-zeta))
	}
}

PROCEDURE rates( v (mV) ) {
	minf = 1 / ( 1 + exp(-(v+cv+vshift)/ck) )
	taum = (1e3) * taumfkt(v)/qt/kt
}

FUNCTION taumfkt( v (mV) ) (s) {
	UNITSOFF

    taumfkt = (0.0002 + 0.0007031 * exp(-((v+30+vshift)/14)^2))				:Raman data
:     taumfkt = (0.00002 + 0.00065 * exp(-((v+vshift)/40)^2))								:data from Biophysical Journal 108,2015: 578-584 David Naranjo
	UNITSON
}

FUNCTION kelvinfkt( t (degC) )  (kelvin) {
	UNITSOFF
	kelvinfkt = 273.19 + t
	UNITSON
}