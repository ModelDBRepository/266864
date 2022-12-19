: Recrd the time at which icap has a peak.
: By Michael Hines 09-2004

NEURON {
 SUFFIX capmax  
 POINTER icap 
 RANGE i,tp,debug,i1,i2,i3,lockit,i1,i2,i3,up
}

PARAMETER { 
	 debug = 0
}

ASSIGNED { 
	 i (milliamp/cm2)  
	 tp  (ms)
	 icap (milliamp/cm2) 
	 i1 
	 i2
	 i3
}
STATE {
      up
      lockit
}

INITIAL { i = 0
	  tp =0
	  lockit = -1
	  up = 0
	  i1= 0
	  i2 = 0
	  i3 = 0
	 }

BREAKPOINT { 
	   SOLVE mx 
}
 
PROCEDURE absmx() {
	   if (fabs(icap) > i) { 
	      i = fabs(icap) 
	      tp  = t
	   }
}

PROCEDURE mxpos() {
	   if (icap > i) { 
	      i = icap 
	      tp  = t
	   }
}


PROCEDURE mx() {
	  if (up == 1 ) { 
	    if (lockit < 0) { 
	   	   if (icap > i) {	
	      	      i = icap 
	      	      tp  =  t
	   	   }
	   	   if (debug == 1 ) {  
	             VERBATIM   	 
	               fprintf(stdout,"Lockit %f\tup: %f\tt: %f\ti: %f\n",lockit,up,t,i);
	             ENDVERBATIM
	          }
            }

          } 
	  if ( (t>10) && (icap>i3) && (i3>i2) && (i2>i1) ) { 
	      up = 1
	  } else { 
	     up =0 
	  }	     	      
	   

	  if ( (icap > 0 ) && ( icap < i*0.75) && ( (t-tp) < 2) ) {  
	       lockit = 1
	  }	       
	  i1 = i2
	  i2 = i3
	  i3 = icap
}




