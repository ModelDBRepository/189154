# $Id: labels.py,v 1.4 2012/09/24 17:35:22 samn Exp $ 
from neuron import h
h.load_file("labels.hoc") # has variables needed by network
hvars=[a for a in dir(h) if h.name_declared(a,2)==5] # 5 if a scalar or double variable
    
CTYP = [s.s for s in h.CTYP]
CTYPi = len(CTYP)
for s in CTYP: globals()[s]=int(h.__getattribute__(s))
E = int(h.E2)
I = int(h.I2)

STYP = [s.s for s in h.STYP]
STYPi = len(STYP)
for s in STYP: globals()[s]=int(h.__getattribute__(s))

try:  # most of these are defined in labels.hoc and syncode.hoc
  # ice - return True iff arg1 is an inhibitory cell
  ice = h.ice
  IsLTS = h.IsLTS # returns true iff arg1 is LTS cell
  SOMA = h.SOMA
  DEND = h.DEND
  AXON = h.AXON
except: print 'something unhandled in labels.py additions 1'
