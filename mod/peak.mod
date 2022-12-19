TITLE peak.mod

COMMENT 
pk: record peak time and peak value of membrane potential
Michael Hausser & Arnd Roth                     25.9.1997
Philipp Vetter	modified last			30.12.1998
Mickey 	London	modified last			30.12.1998

vpeak/tpeak are locked to first peak
dvdt2 is for time t - 2*dt
AP measurements		dvdtmax, dvdtmaxII, onset, vonset, halfwidth, vhalf, vrest, dVdr
electrotonic lengths	X, Xsec
impedance mismatch	Zmismatch, Rmismatch, aZmismatch, aRmismatch, f

Copyright (C) 2001  Philipp Vetter, Arnd Roth and Michael Hausser

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
ENDCOMMENT

UNITS {
	(mV) = (millivolt)
	(Mohms) = (megohms)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX pk
	RANGE tpeak, vpeak, vpeakm,tpeakdvdt
	RANGE dvdtpeak, dvdt2peak, onset_ref, onset, vonset, halfwidth, vhalf, vrest, dvdr
	RANGE Zmismatch, Rmismatch, aZmismatch, aRmismatch, f
	RANGE Zback, Rback, aZback, aRback
	RANGE Zfwd, Rfwd, aZfwd, aRfwd
	RANGE Z, R, aZ, aR
	RANGE Xsec, Xfrc, Xlen, Xo, sign
	RANGE dvdt2,dvdt
	RANGE thresh
}

PARAMETER {
	vhalf (mV)
	v (mV)
	dt (ms)
}

ASSIGNED {
	tpeak (ms)
	tpeakdvdt (ms)
	vpeak (mV)
	vpeakm (mV)
	dvdt (ms)
	dvdt2 (mV/ms*ms)
	dvdt2peak (mV/ms*ms)
	dvdtpeak (mV/ms)
	onset_ref (mV/ms)
	onset (ms)
	v1 (mV)
	v2 (mV)
	v3 (mV)
	vonset (mV)
	halfwidth (ms)
	below
	below_old
	upstroke (ms)
	downstroke (ms)
	dvdr     (mv/micron)
	vrest (mV)
	Rmismatch     (1)
	Zmismatch     (1)
	Rback   (Mohm)
	Zback   (Mohm)
	Rfwd    (Mohm)
	Zfwd    (Mohm)
	R     	(Mohm)
	Z     	(Mohm)


	aRmismatch  (1)
	aZmismatch  (1)
	aRback  (Mohm)
	aZback  (Mohm)
	aRfwd   (Mohm)
	aZfwd   (Mohm)
	aR     	(Mohm)
	aZ     	(Mohm)




	f     (0.001/ms)
	Xo 	(1)
	Xlen    (1)
	Xsec	(1)
	Xfrc	(1)
	sign	(1)
	thresh (nA)
}

INITIAL {
	tpeak = 0 (ms)
	tpeakdvdt = 0 (ms)
	vpeak = -100 (mV)
	vpeakm = -100 (mV)
	onset = 0    (ms)
	dvdtpeak  = 0 (mV/ms)
	dvdt2peak  = 0 (mV/ms)
	downstroke = 0 (ms)
	upstroke   = 0 (ms)
	vrest = v
	check()
}

BREAKPOINT {	SOLVE check  }


PROCEDURE check() {
			if (v > vpeak && vpeakm==-100) { tpeak = t
					 		 vpeak = v }

			if (v+4 < vpeak) { vpeakm = 1 }

			v1        = v2
			v2        = v3
			v3        = v
		        dvdt	  = (v3 - v2)/(dt)
			dvdt2	  = (v3 - 2*v2 +v1)/(dt*dt)
			if (t < 1)  {
			dvdt2 =0
			}
			

			if (dvdt > dvdtpeak) 			      	{ 
			   	   					  dvdtpeak = dvdt 
			   	   					  tpeakdvdt = t	     	  
									}
			if (dvdt2 > dvdt2peak) 			      	{ dvdt2peak = dvdt2 }

			if (dvdt > onset_ref && onset == 0 && t > 1) 	{ onset  = t-dt
							                  vonset = v2 }
			below     = 0	
			if (vhalf > v) {below = 1}
			if (below == 0 && below_old == 1 && upstroke == 0)   {upstroke = t}
			if (below == 1 && below_old == 0 && downstroke == 0) {downstroke = t}

			halfwidth = downstroke - upstroke
			below_old = below
}
