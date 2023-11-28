UNITS {	
  (mollar) = (1/liter)
  (M)      = (mollar)
  (mM)     = (millimollar)
  (mA)     = (milliamp)
  (mV)     = (millivolt)
  (mS)     = (millisiemens)
}

NEURON {
  SUFFIX calts
  USEION ca READ ica, cai WRITE cai
  RANGE alpha, tau
}
	
PARAMETER {
  alpha = 0.002 (cm2-M/mA-ms)
  tau   = 80    (ms)
}
    
ASSIGNED { ica (mA/cm2) }

INITIAL { cai = 0 }

STATE { cai (mM) }

BREAKPOINT { SOLVE states METHOD cnexp }

DERIVATIVE states { cai' = -(1000) * alpha * ica - cai/tau }
