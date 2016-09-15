UNITS {
  (mA)     = (milliamp)
  (mV)     = (millivolt)
  (mS)     = (millisiemens)
  (mollar) = (1/liter)
  (mM)     = (millimollar)
}

NEURON {
  SUFFIX kcalts
  USEION k WRITE ik
  USEION ca READ cai
  RANGE gkca,ek,kd
}
	
PARAMETER {
  gkca =  10 (mS/cm2)
  ek   = -90 (mV)
  kd   =  30 (mM)
}
    
ASSIGNED {    
  cai (mM) 
  v   (mV)
  ik  (mA/cm2) 
}

PROCEDURE iassign () { ik  = (1e-3) * gkca * cai/(cai+kd) * (v-ek) }

INITIAL {
  iassign()
}

BREAKPOINT { iassign() }
