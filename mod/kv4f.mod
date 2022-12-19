TITLE Voltage-gated potassium channel from Kv4 subunits

COMMENT

NEURON implementation of a potassium channel from Kv4 subunits
Kv4 activation from Sacco inactivation from SD
Yunliang Zang April 16th 2015
activation from
Channel Density Distributions Explain Spiking Variability in the Globus Pallidus: A Combined Physiology and Computer Simulation Database Approach

ENDCOMMENT

NEURON {
	SUFFIX Kv4
	USEION k READ ek WRITE ik
	RANGE gk, gbar, ik,vshift
:   GLOBAL ninf, taun, hinf, tauh
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
	F = 9.6485e4 (coulombs)
	R = 8.3145 (joule/kelvin)
	can = 0.15743 (1/ms)
	cvan = 57 (mV)
	ckan = -32.19976 (mV)
	cbn = 0.15743 (1/ms)
	cvbn = 57 (mV)
	ckbn = 37.51346 (mV)

	
	cah = 0.01342 (1/ms)
	cvah = 60 (mV)
	ckah = -7.86476 (mV)
	cbh = 0.04477 (1/ms)
	cvbh = 54 (mV)
	ckbh = 11.3615 (mV)
	
	vh = -75.30348 (mV)
	kh = -6.06329 (mV)
	ki = 150 (mM)			:from Stephane
	ko = 2.5 (mM)   
}

PARAMETER {
	v (mV)
	celsius (degC)
	vshift = 0
	gbar = 0.0039 (mho/cm2)   <0,1e9>
}

ASSIGNED {
	ik (mA/cm2) 
	ek (mV)
	gk (mho/cm2)
	g (coulombs/cm3)
	
	T (kelvin)
	qt
	E (volt)
	zeta

	ninf
	taun (ms)
	alphan (1/ms)
	betan (1/ms)
	alphah (1/ms)
	betah (1/ms)	

	hinf
:	h1inf
:	h2inf
	tauh (ms)
:	tauh2 (ms)    
}

STATE { n h }

INITIAL {
	T = kelvinfkt (celsius)
	qt = q10^((celsius-23 (degC))/10 (degC))
	rates(v)
	n = ninf
	h = hinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
    gk = gbar * n*n*n*n*h
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
	
: activation from Jager	
	ninf = 1.0 / (1.0 + exp((-49 - v)/12.5))

	taun = 1/((alphan+betan)*qt)
	alphah = alphahfkt(v)
	betah = betahfkt(v)

	hinf = 1/(1+exp((v-(vh-vshift))/-kh))
	tauh =20/qt
	g = ghk(v, ki, ko, 1)
}

FUNCTION ghk( v (mV), ki (mM), ko (mM), z )  (coulombs/cm3) {
    E = (1e-3) * v
      zeta = (z*F*E)/(R*T)


    if ( fabs(1-exp(-zeta)) < 1e-6 ) {
        ghk = (1e-6) * (z*F) * (ki - ko*exp(-zeta)) * (1 + zeta/2)
    } else {
        ghk = (1e-6) * (z*zeta*F) * (ki - ko*exp(-zeta)) / (1-exp(-zeta))
    }
}


FUNCTION alphanfkt(v (mV)) (1/ms) {
	alphanfkt = can * exp(-(v+cvan)/ckan) 
}

FUNCTION betanfkt(v (mV)) (1/ms) {
	betanfkt = cbn * exp(-(v+cvbn)/ckbn)
}

FUNCTION kelvinfkt( t (degC) )  (kelvin) {
    kelvinfkt = 273.19 + t
}
FUNCTION alphahfkt(v (mV))  (1/ms) {
	alphahfkt = cah / (1+exp(-(v+cvah)/ckah))
}

FUNCTION betahfkt(v (mV))  (1/ms)  {
	betahfkt = cbh / (1+exp(-(v+cvbh)/ckbh))
}
