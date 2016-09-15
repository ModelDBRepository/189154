from neuron import h

# metabotropic glutamate receptor
class SynapsemGLUR:
  def __init__(self,sect,loc):
    self.syn = h.mGLUR(loc, sec=sect)

# AMPA synapse with calcium influx -- mechanism defined in mod/ampa_forti.mod
class SynapseAMPACA:
  def __init__(self, sect, loc, e):
    self.syn		= h.AmpaSyn(loc, sec=sect)
    self.syn.e		= e 

# NMDA synapse with calcium influx -- mechanism defined in mod/nmda_andr.mod
class SynapseNMDACA:
  def __init__(self, sect, loc, e):
    self.syn		= h.NmdaSyn(loc, sec=sect)
    self.syn.e		= e 

class Synapse:
  def __init__(self, sect, loc, tau1, tau2, e):
    self.syn		= h.MyExp2SynBB(loc, sec=sect)
    self.syn.tau1	= tau1
    self.syn.tau2	= tau2
    self.syn.e		= e 
		
class SynapseNMDA:
  def __init__(self, sect, loc, tau1NMDA, tau2NMDA, r, e):
    self.syn			= h.MyExp2SynNMDABB(loc, sec=sect)
    self.syn.tau1NMDA	= tau1NMDA
    self.syn.tau2NMDA	= tau2NMDA 
    self.syn.r			= r
    self.syn.e			= e 

# gabab based on 1995 PNAS paper by Destexhe
class SynapseGABAB:
  def __init__(self, sect, loc):
    self.syn = h.GABAB(loc, sec=sect)

class SynapseSTDP:
  def __init__(self, sect, loc, tau, e, dtau, ptau, d, p):
    self.syn	= h.ExpSynSTDP(loc, sec=sect)
    self.syn.tau    = tau
    self.syn.e     	= e 
    self.syn.dtau	= dtau
    self.syn.ptau	= ptau
    self.syn.d      = d
    self.syn.p      = p
