: Calcium ion accumulation with radial and longitudinal diffusion and pump

NEURON {
  SUFFIX cdpAIS
  USEION ca READ cao, cai, ica WRITE cai
  RANGE ica_pmp,ca1,ca2,ca3,ca4,ca5,ca9
:RANGE pump_0  
GLOBAL vrat, TotalPump
    : vrat must be GLOBAL--see INITIAL block
    : however TotalBuffer and TotalPump may be RANGE
:    THREADSAFE
}

DEFINE Nannuli 10

UNITS {
	(mol)   = (1)
	(molar) = (1/liter)
	(mM)    = (millimolar)
	(um)    = (micron)
	(mA)    = (milliamp)
	FARADAY = (faraday)  (10000 coulomb)
	PI      = (pi)       (1)
}

PARAMETER {
    v
	celsius =34     (degC)
        
	:cainull =2.5e-4 (mM)
	cainull = 45e-6 (mM)
        mginull =.59    (mM)

        DCa     = .233  (um2/ms)
	Dbtc 	= 0.007 (um2/ms)
       Ddmnpe = 0.08	(um2/ms)
	
	Dcbd1   = .028  (um2/ms)
        Dcbd2   = 0     (um2/ms)
        Dpar    = .043  (um2/ms)

:	values for benzothiazole coumarin (BTC)
:	BTCnull = 0	(mM)
:	b1 = 5.33	(/ms mM)
:	b2 = 0.08	(/ms)

:	values for caged compound DMNPE-4
:	DMNPEnull = 0	(mM)
:	c1 = 5.63	(/ms mM)
:	c2 = 0.107e-3	(/ms)

:       values for Calbindin (2 high and 2 low affinity binding sites)

:        CBnull=	.16             (mM)
         CBnull=	.08             (mM)       
        nf1   =43.5           (/ms mM)
        nf2   =3.58e-2        (/ms)
        ns1   =5.5            (/ms mM)
        ns2   =0.26e-2        (/ms)

:       values for Parvalbumin

:        PVnull  = 0.08           (mM)
        PVnull  = 0.04           (mM)
        m1    = 1.07e2        (/ms mM)
        m2    = 9.5e-4                (/ms)
        p1    = 0.8           (/ms mM)
        p2    = 2.5e-2                (/ms)        
    
  	kpmp1    = 3e3       (/mM-ms)
  	kpmp2    = 1.75e1   (/ms)
  	kpmp3    = 7.255e1  (/ms)
  : to eliminate pump, set TotalPump to 0 in hoc
	TotalPump = 1e-15
:	TotalPump = 1e-8
	beta  = 1(1)           :introducing beta to take care of other ER mechanisms(SERCA and leak channel density)
    vmax =0.1
    :Kp = 0.27e-3 (mM)
    Kp = 2.7e-3 (mM)
}

ASSIGNED {
	diam      (um)
	ica       (mA/cm2)
	ica_pmp   (mA/cm2)
	ibg     :background calcium current
:	ica_pmp_last   (mA/cm2)
	parea     (um)     : pump area per unit length
	cai       (mM)
	ca1
	ca2
	ca3
	ca4
	ca5
	ca9
	mgi	(mM)	
	vrat[Nannuli]  (1) : dimensionless
                     : numeric value of vrat[i] equals the volume 
                     : of annulus i of a 1um diameter cylinder
                     : multiply by diam^2 to get volume per um length
	
}

CONSTANT { cao = 2	(mM) }

STATE {
	: ca[0] is equivalent to cai
	: ca[] are very small, so specify absolute tolerance
	: let it be ~1.5 - 2 orders of magnitude smaller than baseline level
	ca[Nannuli]		(mM)
	mg[Nannuli]		(mM)	<1e-7>

        CB[Nannuli]		(mM)
        CB_f_ca[Nannuli]	(mM)
        CB_ca_s[Nannuli]	(mM)
        CB_ca_ca[Nannuli]	(mM)

        iCB[Nannuli]		(mM)
        iCB_f_ca[Nannuli]	(mM)
        iCB_ca_s[Nannuli]	(mM)
        iCB_ca_ca[Nannuli]	(mM)

        PV[Nannuli]		(mM)
        PV_ca[Nannuli]		(mM)
        PV_mg[Nannuli]		(mM)
	
	pump			(mol/cm2) <1e-15>
	pumpca			(mol/cm2) <1e-15>
}

BREAKPOINT {
	SOLVE state METHOD sparse
:	ica_pmp_last = ica_pmp
:	ica = ica_pmp
}

LOCAL factors_done

INITIAL {
	if (factors_done == 0) {  : flag becomes 1 in the first segment
		factors_done = 1       :   all subsequent segments will have
		factors()              :   vrat = 0 unless vrat is GLOBAL
	}
	FROM i=0 TO Nannuli-1 {
		ca[i] = cainull
		mg[i] = mginull


		CB[i] = 0.8*ssCB( kdf(), kds())   
	        CB_f_ca[i] = 0.8*ssCBfast( kdf(), kds())
       	 	CB_ca_s[i] = 0.8*ssCBslow( kdf(), kds())
        	CB_ca_ca[i] = 0.8*ssCBca( kdf(), kds())

        	iCB[i] = 0.2*ssCB( kdf(), kds())
        	iCB_f_ca[i] = 0.2*ssCBfast( kdf(), kds())
        	iCB_ca_s[i] = 0.2*ssCBslow( kdf(), kds())
        	iCB_ca_ca[i] = 0.2*ssCBca(kdf(), kds())

        	PV[i] = ssPV( kdc(), kdm())
        	PV_ca[i] = ssPVca(kdc(), kdm())
        	PV_mg[i] = ssPVmg(kdc(), kdm())

		
	}
  	parea = PI*diam
	ica = 0
	ica_pmp = 0
:	ica_pmp_last = 0
	pump = TotalPump
	pumpca = 0
}

LOCAL frat[Nannuli]  : scales the rate constants for model geometry

PROCEDURE factors() {
	LOCAL r, dr2
	r = 1/2                : starts at edge (half diam)
	dr2 = r/(Nannuli-1)/2  : full thickness of outermost annulus,
                         : half thickness of all other annuli
	vrat[0] = 0
	frat[0] = 2*r
	FROM i=0 TO Nannuli-2 {
		vrat[i] = vrat[i] + PI*(r-dr2/2)*2*dr2  : interior half
		r = r - dr2
		frat[i+1] = 2*PI*r/(2*dr2)  : outer radius of annulus
                                : div by distance between centers
		r = r - dr2
    		vrat[i+1] = PI*(r+dr2/2)*2*dr2  : outer half of annulus
  	}
}

LOCAL dsq, dsqvol  : can't define local variable in KINETIC block
                   :   or use in COMPARTMENT statement

KINETIC state {
  COMPARTMENT i, diam*diam*vrat[i] {ca mg CB CB_f_ca CB_ca_s CB_ca_ca iCB iCB_f_ca iCB_ca_s iCB_ca_ca PV PV_ca PV_mg}
  COMPARTMENT (1e10)*parea {pump pumpca}
	:pump
:	~ ca[0] + pump <-> pumpca  (kpmp1*parea*(1e10), kpmp2*parea*(1e10))
:	~ pumpca <-> pump   (kpmp3*parea*(1e10), 0)
:  	CONSERVE pump + pumpca = TotalPump * parea * (1e10)

:	ica_pmp = 2*FARADAY*(f_flux - b_flux)/parea

	: all currents except pump
	: ica is Ca efflux
	~ ca[0] << (-ica*PI*diam/(2*FARADAY))

	:RADIAL DIFFUSION OF ca, mg and mobile buffers

    FROM i=0 TO Nannuli-1 {
    ~ ca[i] << (-beta*vmax*vrat[i]*ca[i] / (ca[i] + kpmp2/kpmp1))

    }
    
	FROM i=0 TO Nannuli-2 {
		~ ca[i] <-> ca[i+1]	(DCa*frat[i+1], DCa*frat[i+1])
		~ mg[i] <-> mg[i+1]	(DCa*frat[i+1], DCa*frat[i+1])
		~ CB[i] <-> CB[i+1]	(Dcbd1*frat[i+1], Dcbd1*frat[i+1])
		~ CB_f_ca[i] <-> CB_f_ca[i+1]	(Dcbd1*frat[i+1], Dcbd1*frat[i+1])
		~ CB_ca_s[i] <-> CB_ca_s[i+1]	(Dcbd1*frat[i+1], Dcbd1*frat[i+1])
		~ CB_ca_ca[i] <-> CB_ca_ca[i+1]	(Dcbd1*frat[i+1], Dcbd1*frat[i+1])
		~ PV[i] <-> PV[i+1]	(Dpar*frat[i+1], Dpar*frat[i+1])
		~ PV_ca[i] <-> PV_ca[i+1]	(Dpar*frat[i+1], Dpar*frat[i+1])
		~ PV_mg[i] <-> PV_mg[i+1] 	(Dpar*frat[i+1], Dpar*frat[i+1])
	}
	dsq = diam*diam
	FROM i=0 TO Nannuli-1 {
		dsqvol = dsq*vrat[i]
		:Calbindin	
		~ ca[i] + CB[i] <-> CB_ca_s[i] (nf1*dsqvol, nf2*dsqvol)
	       	~ ca[i] + CB[i] <-> CB_f_ca[i] (ns1*dsqvol, ns2*dsqvol)
        	~ ca[i] + CB_f_ca[i] <-> CB_ca_ca[i] (nf1*dsqvol, nf2*dsqvol)
        	~ ca[i] + CB_ca_s[i] <-> CB_ca_ca[i] (ns1*dsqvol, ns2*dsqvol)

        	~ ca[i] + iCB[i] <-> iCB_ca_s[i] (nf1*dsqvol, nf2*dsqvol)
        	~ ca[i] + iCB[i] <-> iCB_f_ca[i] (ns1*dsqvol, ns2*dsqvol)
        	~ ca[i] + iCB_f_ca[i] <-> iCB_ca_ca[i] (nf1*dsqvol, nf2*dsqvol)
        	~ ca[i] + iCB_ca_s[i] <-> iCB_ca_ca[i] (ns1*dsqvol, ns2*dsqvol)


		:Paravalbumin
        	~ ca[i] + PV[i] <-> PV_ca[i] (m1*dsqvol, m2*dsqvol)
        	~ mg[i] + PV[i] <-> PV_mg[i] (p1*dsqvol, p2*dsqvol)
	}



  	cai = ca[0]
  	ca1 = ca[1]
  	ca2 = ca[2]
  	ca3 = ca[3]
  	ca4 = ca[4]
  	ca5 = ca[5]
  	ca9 = ca[9]
:  	if (ca[0]<cainull){: keep minimum
:	   cai=cainull }
	mgi = mg[0]
}


FUNCTION ssCB( kdf(), kds()) (mM) {
	ssCB = CBnull/(1+kdf()+kds()+(kdf()*kds()))
}
FUNCTION ssCBfast( kdf(), kds()) (mM) {
	ssCBfast = (CBnull*kds())/(1+kdf()+kds()+(kdf()*kds()))
}
FUNCTION ssCBslow( kdf(), kds()) (mM) {
	ssCBslow = (CBnull*kdf())/(1+kdf()+kds()+(kdf()*kds()))
}
FUNCTION ssCBca(kdf(), kds()) (mM) {
	ssCBca = (CBnull*kdf()*kds())/(1+kdf()+kds()+(kdf()*kds()))
}
FUNCTION kdf() (1) {
	kdf = (cainull*nf1)/nf2
}
FUNCTION kds() (1) {
	kds = (cainull*ns1)/ns2
}
FUNCTION kdc() (1) {
	kdc = (cainull*m1)/m2
}
FUNCTION kdm() (1) {
	kdm = (mginull*p1)/p2
}
FUNCTION ssPV( kdc(), kdm()) (mM) {
	ssPV = PVnull/(1+kdc()+kdm())
}
FUNCTION ssPVca( kdc(), kdm()) (mM) {
	ssPVca = (PVnull*kdc)/(1+kdc+kdm)
}
FUNCTION ssPVmg( kdc(), kdm()) (mM) {
	ssPVmg = (PVnull*kdm())/(1+kdc()+kdm())
}


FUNCTION u (x, th) {
  	if (x<th) {
    		u = 1
  	} else {
    		u = 0
  	}
}
