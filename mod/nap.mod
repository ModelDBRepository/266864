: nap.mod is a persistent Na+ current from
: Baker 2005, parameter assignments and formula's from page 854

NEURON {
	SUFFIX nap
    USEION na READ ena WRITE ina
	RANGE gbar,ina
:    THREADSAFE
}

UNITS {
	(S) = (siemens)
	(mV) = (millivolts)
	(mA) = (milliamp)
}
CONSTANT {
q10  =2.7
}
PARAMETER {
	gbar = 2.2630e-04 :3.7(nS)/1635(um^2)
:	ena= 65 (mV)

	A_amp = 17.235 (/ms) : A for alpha m persis
	B_amp = 27.58 (mV)
	C_amp = -11.47 (mV)

	A_bmp = 17.235 (/ms) : A for beta m persis
	B_bmp = 86.2 (mV)
	C_bmp = 19.8 (mV)
}

ASSIGNED {
	v	(mV) : NEURON provides this
	i	(mA/cm2)
	g	(S/cm2)
	tau_m	(ms)
	minf
	hinf
	ena
	ina
	qt
}

STATE { m h }

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar * m^3
	ina = g * (v-ena)
}

INITIAL {
qt = q10^((celsius-22 (degC))/10 (degC))
	: assume that equilibrium has been reached
	m = alpham(v)/(alpham(v)+betam(v))
}

DERIVATIVE states {
	rates(v)
	m' = (minf - m)/tau_m
}

FUNCTION alpham(Vm (mV)) (/ms) {
	alpham=A_amp/(1+exp((Vm+B_amp)/C_amp))
}

FUNCTION betam(Vm (mV)) (/ms) {
	betam=A_bmp/(1+exp((Vm+B_bmp)/C_bmp))
}

FUNCTION rates(Vm (mV)) (/ms) {
	tau_m = 1.0 / (alpham(Vm) + betam(Vm))/qt
	minf = 1.0/(1+exp(-(Vm+66)/5))
}












