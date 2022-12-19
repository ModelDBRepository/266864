TITLE ...just to store peak membrane voltage
: M.Migliore June 2001
: T Morse February 2010 added times of occurrence

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	v (mV)
}


NEURON {
	SUFFIX ds
        RANGE vmax, tmax,vamp
}

ASSIGNED {
	vmax
	tmax
	vmin
	vamp
}

INITIAL {
:	vmax=v
    vmax = -90
    vmin = -30
    vamp = 0
}


BREAKPOINT {
    if (t>50) {
    if (v<vmin) {vmin=v}
	if (v>vmax) {vmax=v tmax=t}
	vamp = vmax-vmin
	}
	else {
	vmax = -70
	tmax=0
	}
	
}
