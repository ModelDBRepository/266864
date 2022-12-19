TITLE Low threshold calcium current Cerebellum Purkinje Cell Model

COMMENT

Q10 is estimated from this work, Temperature dependence of T-type Calcium channel gating, NEUROSCIENCE
written by Yunliang Zang according to the data provided by Stephane Diudone, compared with the summarised data from stephane,
T type calcium channels has two gates. so the activation curve was refitted.
The junction potential is -6.6 mV
It does not work even changing it back to cai
April 16th, 2015
This version does not contribute to the calcium concentration and BK together with SK. 
ENDCOMMENT


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
        SUFFIX CaT3_1
:        USEION ca READ cai, cao WRITE ica VALENCE 2
:        NONSPECIFIC_CURRENT i
	USEION ca READ cai,cao
	USEION Ca WRITE iCa VALENCE 2
        RANGE g, pcabar, minf, taum, hinf, tauh
    	RANGE iCa, m ,h
:    THREADSAFE
    }

UNITS {
        (molar) = (1/liter)
        (mV) =  (millivolt)
        (mA) =  (milliamp)
        (mM) =  (millimolar)

}

CONSTANT {
    q10 = 1.0913        :estimate from Iftinca
	F = 9.6485e4 (coulombs)
	R = 8.3145 (joule/kelvin)
}

PARAMETER {
        v               (mV)
        celsius (degC)
        eca (mV)
	pcabar  = 2.5e-4 (cm/s)
        cai = 1e-4  (mM)           : adjusted for eca=120 mV
	cao = 2  (mM)
	
	v0_m_inf = -42.206 (mV)
	v0_h_inf = -75.118 (mV)
	vshift = -6.6			:liquid junction potential

	k_m_inf = -4.7056 (mV)
	k_h_inf = 6.4635  (mV)
	
	C_tau_m = 1.2757
	A_tau_m = -2.3199
	B_tau_m = 2.5712
	v0_tau_m1 = -48.048 (mV)
	v0_tau_m2 = -28.386 (mV)
	k_tau_m1 = 30.655 (mV)
	k_tau_m2 = 9.6306 (mV)
	
	C_tau_h = 0.0076
	A_tau_h = 0.17746
	B_tau_h = 0.13402
	v0_tau_h1 = -58.535 (mV)
	v0_tau_h2=-101.436
	k_tau_h1 = 6.2692 (mV)
	k_tau_h2 = -5.5845 (mV)
}
    

STATE {
        m h
}

ASSIGNED {
        iCa     (mA/cm2)
	g        (coulombs/cm3) 
        minf
        taum   (ms)
        hinf
        tauh   (ms)
        qt
	T (kelvin)
	E (volt)
	zeta
}

BREAKPOINT {
	SOLVE castate METHOD cnexp 

       iCa = (1e3) *pcabar*m*m *m*h * g

}

DERIVATIVE castate {
        evaluate_fct(v)

        m' = (minf - m) / taum
        h' = (hinf - h) / tauh
}

FUNCTION ghk2( v (mV), ci (mM), co (mM), z )  (coulombs/cm3) {
    E = (1e-3) * v
      zeta = (z*F*E)/(R*T)


    if ( fabs(1-exp(-zeta)) < 1e-6 ) {
        ghk2 = (1e-6) * (z*F) * (ci - co*exp(-zeta)) * (1 + zeta/2)
    } else {
        ghk2 = (1e-6) * (z*zeta*F) * (ci - co*exp(-zeta)) / (1-exp(-zeta))
    }
}


UNITSOFF
INITIAL {
	
	T = kelvinfkt (celsius)
	    qt = q10^((celsius-32 (degC))/10 (degC))
        evaluate_fct(v)
        m = minf
        h = hinf
}

PROCEDURE evaluate_fct(v(mV)) { 

        minf = 1.0 / ( 1 + exp((v  - v0_m_inf-vshift)/k_m_inf) )^(1/3)
        
        hinf = 1.0 / ( 1 + exp((v - v0_h_inf-vshift)/k_h_inf) )

	taum = 1/( C_tau_m + A_tau_m / (1+exp((v0_tau_m1-v-vshift)/ k_tau_m1))+ B_tau_m/ (1+exp((v0_tau_m2-v-vshift)/k_tau_m2)))/qt

	tauh = 1/( C_tau_h + A_tau_h / (1+exp((v0_tau_h1-v-vshift)/ k_tau_h1))+ B_tau_h/ (1+exp((v0_tau_h2-v-vshift)/k_tau_h2)))/qt

	g = ghk2(v-vshift, cai, cao, 2)
}

FUNCTION kelvinfkt( t (degC) )  (kelvin) {
    kelvinfkt = 273.19 + t
}

UNITSON
