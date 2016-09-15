import sys
from pyinit import *
from labels import *
from math import exp
h.celsius = 37
h.load_file("pywrap.hoc")
from conf import *
# determine config file name
def setfcfg ():
  fcfg = "netcfg.cfg" # default config file name
  for i in xrange(len(sys.argv)):
    if sys.argv[i].endswith(".cfg") and os.path.exists(sys.argv[i]):
      fcfg = sys.argv[i]
  return fcfg

fcfg=setfcfg() # config file name
dconf = readconf(fcfg)
taurcada = dconf['taurcada']
h.cac_hcnwino = 0.006
h.k4_hcnwino = dconf['iark4'] 
ihginc = h.ginc_hcnwino = dconf['ihginc'];
recdt = dconf['recdt']
recvdt = dconf['recvdt']
erevh = dconf['erevh']
spaceum = dconf['spaceum']
h_lambda = dconf['h_lambda']
h_gbar = dconf['h_gbar'] # for E cells
fs_h_gbar = dconf['fs_h_gbar'] # 
lts_h_gbar = dconf['lts_h_gbar'] # 
cagk_gbar = dconf['cagk_gbar'] # 
ikc_gkbar = dconf['ikc_gkbar'] # 
cabar = dconf['cabar'] # used for E cells
lts_cabar = dconf['lts_cabar']
tau1NMDAEE=15; tau2NMDAEE=150;
tau1NMDAEI=15; tau2NMDAEI=150;
nax_gbar = dconf['nax_gbar'] 
kdr_gbar = dconf['kdr_gbar'] 
kap_gbar = dconf['kap_gbar'] 
kdmc_gbar = dconf['kdmc_gbar'] 
km_gmax = dconf['km_gmax'] 
##

from syn import *

# if rdt > 0 use fixed interval for recording, else let cvode determine it
def saferecord (var, rdt):
  if rdt > 0.0:
    vrec = h.Vector(h.tstop/rdt + 1)
    vrec.record(var,rdt)
  else:
    vrec = h.Vector()
    vrec.record(var)
  return vrec
		
###############################################################################
# General Cell
###############################################################################
class Cell:
  "General cell"
  def __init__ (self,x,y,z,ID,ty):
    self.x=x
    self.y=y
    self.z=z
    self.ID=ID
    self.ty = ty
    self.snames = [] # list of section names
    self.all_sec = []
    self.add_comp('soma',True)
    self.set_morphology()
    self.set_conductances()
    self.set_synapses()
    self.set_inj()
    
  # get number of outgoing connections
  def set_morphology (self): pass			
  def set_conductances (self): pass
  def set_synapses (self): pass
  def set_inj (self): self.somaInj = h.IClamp(0.5, sec=self.soma)	
		
  def add_comp (self, name, rec):
    self.snames.append( name )
    self.__dict__[name] = h.Section()
    self.all_sec.append(self.__dict__[name])
    if rec:    # Record voltage
      self.__dict__[name+"_volt"] = saferecord(self.__dict__[name](0.5)._ref_v, recvdt)
      self.__dict__[name+"_volt"].label(name+"_volt")
      
###############################################################################
# Soma-targeting interneuron (fast-spiking Basket Cell -- Bas)
###############################################################################
class Bas (Cell):
  "Basket cell"	
  def set_morphology(self):
    total_area = 10000 # um2
    self.soma.nseg  = 1
    self.soma.cm    = 1      # uF/cm2
    diam = sqrt(total_area) # um
    L    = diam/pi  # um			
    h.pt3dclear(sec=self.soma)
    h.pt3dadd(self.x, self.y, self.z,   diam, sec=self.soma)
    h.pt3dadd(self.x, self.y, self.z+L, diam, sec=self.soma)
			
  def set_conductances(self): # Bas
    cap         = 1.0
    rall        = 150.0
    rm          = 10e3 
    Vrest       = -79.8
    p_ek          = -85.0 
    p_ena        = 55.0 
    sh_nax = 0.0
    sec = self.soma
    sec.insert('k_ion')
    sec.insert('na_ion')
    #sec.insert('ca_ion')
    sec.insert('pas') # passive     
    sec.insert('nax') # Na current      
    sec.insert('kdr') # K delayed rectifier current
    # erev
    sec.ek = p_ek # K+ current reversal potential (mV)
    sec.ena = p_ena # Na+ current reversal potential (mV)
    # passive
    sec.g_pas = 1.0/rm
    sec.Ra = rall
    sec.cm = cap
    sec.e_pas = Vrest
    # Na
    sec.gbar_nax = nax_gbar
    sec.sh_nax = sh_nax
    # KDR
    sec.gbar_kdr = kdr_gbar
    self.soma.insert('HCN1')
    self.soma(0.5).HCN1.gbar = fs_h_gbar
	   
  def set_synapses(self):
    self.somaGABAf=Synapse(sect=self.soma,loc=0.5,tau1=0.07,tau2=9.1,e=-80);#self.dSy['somaGABAf']=self.somaGABAf;
    self.somaGABAss=Synapse(sect=self.soma,loc=0.5,tau1=20,tau2=40,e=-80);#self.dSy['somaGABAss']=self.somaGABAss;
    self.somaAMPA=Synapse(sect=self.soma,loc=0.5,tau1=0.05,tau2=5.3,e=0);#self.dSy['somaAMPAf']=self.somaAMPAf;
    self.somaNMDA=SynapseNMDA(sect=self.soma,loc=0.5, tau1NMDA=tau1NMDAEI,tau2NMDA=tau2NMDAEI,r=1,e=0);
		
###############################################################################
# Dendrite-targeting interneuron (LTS Cell)
###############################################################################
class Lts (Cell):
  "LTS cell"   
  def set_morphology(self):
    total_area = 10000 # um2
    self.soma.nseg  = 1
    self.soma.cm    = 1      # uF/cm2
    diam = sqrt(total_area) # um
    L    = diam/pi  # um
    h.pt3dclear(sec=self.soma)
    h.pt3dadd(self.x, self.y, self.z,   diam, sec=self.soma)
    h.pt3dadd(self.x, self.y, self.z+L, diam, sec=self.soma)
	
  def set_conductances(self): # LTS
    cap         = 1.0
    rall        = 150.0
    rm          = 10e3 
    Vrest       = -79.8
    p_ek          = -85.0 
    p_ena        = 55.0 
    sh_nax = 0.0
    sec = self.soma
    sec.insert('k_ion')
    sec.insert('na_ion')
    sec.insert('pas') # passive     
    sec.insert('nax') # Na current      
    sec.insert('kdr') # K delayed rectifier current
    # erev
    sec.ek = p_ek # K+ current reversal potential (mV)
    sec.ena = p_ena # Na+ current reversal potential (mV)
    # passive
    sec.g_pas = 1.0/rm
    sec.Ra = rall
    sec.cm = cap
    sec.e_pas = Vrest
    # Na
    sec.gbar_nax = nax_gbar
    sec.sh_nax = sh_nax
    # KDR
    sec.gbar_kdr = kdr_gbar
    # ca-related 
    sec.insert('icalts')
    sec(0.5).icalts.gca = lts_cabar
    sec.insert('kcalts')
    sec.insert('ihlts')
    sec(0.5).ihlts.gh = lts_h_gbar
    sec.insert('calts') # calcium extrusion
    sec(0.5).calts.tau = taurcada
    
  def set_synapses(self):
    self.somaGABAf 	= Synapse(sect=self.soma, loc=0.5, tau1=0.07, tau2=9.1, e=-80)
    self.somaGABAss	= Synapse(    sect=self.soma, loc=0.5, tau1=20,	  tau2=40, e=-80)#originally for septal input
    self.somaAMPA 	= Synapse(    sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, e=0)
    self.somaNMDA 	= SynapseNMDA(sect=self.soma, loc=0.5, tau1NMDA=tau1NMDAEI, tau2NMDA=tau2NMDAEI, r=1, e=0)

LTS = Lts
FS = Bas
		
###############################################################################
# Pyramidal Cell
###############################################################################
class PyrAdr (Cell):
  "Pyramidal cell"
  def __init__(self,x,y,z,ID,ty):
    Cell.__init__(self,x,y,z,ID,ty)
    self.set_props()
    lrec = ['soma','Adend3']

  def set_morphology(self):
    self.add_comp('Bdend',True)
    self.add_comp('Adend1',False)
    self.add_comp('Adend2',False)
    self.add_comp('Adend3',True)
    self.apic = [self.Adend1, self.Adend2, self.Adend3]
    self.basal = [self.Bdend]
    sec = self.soma; sec.L = 20.0; sec.diam = 20.0
    if self.ty == E5R or self.ty == E5B or self.ty == E5P: apicL = 300.0
    else: apicL = 150.0
    #else: apicL = 300.0
    for sec in self.apic:
      sec.L = apicL; sec.diam = 2.0
    self.Bdend.L = 200.0; self.Bdend.diam = 2.0

    self.Bdend.connect(self.soma,    0, 0)
    self.Adend1.connect(self.soma,   1, 0)
    self.Adend2.connect(self.Adend1, 1, 0)
    self.Adend3.connect(self.Adend2, 1, 0)

    if spaceum > 0.0:
      for sec in self.all_sec:
        ns = int(sec.L / spaceum)
        if ns % 2 == 0: ns += 1
        sec.nseg = ns

  def set_props (self): # PYR
    Vrest       = -79.8 
    h.v_init = -60.0
    #h.v_init = -79.8 # -70 # -75 # -79.8 # Vrest # -79.8
    # passive properties
    cap         = 1.0
    rall        = 150.0
    rm          = 10e3 
    # Na, K reversal potentials calculated from
    # internal and external solutions via Nernst equation
    p_ek          = -85.0 
    p_ena        = 55.0 
    # h-current 
    #h.erev_h      = -42.0
    gbar_h      = h_gbar 
    # d-current 
    kdmc_gbar_somam = 20
    # na,k 
    sh_nax = 0.0
    gbar_nax    = nax_gbar
    nax_gbar_somam = 5
    kdr_gbar_somam = 5
    # A few kinetic params changed vis-a-vis kdr_BS.mod defaults:
    h.a0n_kdr     = 0.0075 # def 0.02
    h.nmax_kdr    = 20.0 # def 2
    sh_kap = 0.0
    kap_gbar_somam = 5
    # A few kinetic params changed from kap_BS.mod defaults:
    h.vhalfn_kap  = 35.0 # def 11
    h.nmin_kap    = 0.4 # def 0.1
    h.lmin_kap    = 5.0 # def 2
    h.tq_kap      = -45.0 # def -40
    # other ion channel parameters
    cal_gcalbar = cabar 
    can_gcanbar = cabar 
    cat_gcatbar = cabar 
    calginc = 1.0 # 2.0 - middle might need to get more but can leave out
    cal_gbar_somam = can_gbar_somam = cat_gbar_somam = 0.1
    cal_gbar_bdendm = can_gbar_bdendm = cat_gbar_bdendm = 0.25
    ikc_gbar_dendm = 0.25
    for sec in self.all_sec:
      # erev
      sec.ek = p_ek # K+ current reversal potential (mV)
      sec.ena = p_ena # Na+ current reversal potential (mV)
      # passive
      sec.g_pas = 1.0/rm
      sec.Ra = rall
      sec.cm = cap
      sec.e_pas = Vrest
      # Ih
      sec.ehwino = erevh
      for seg in sec:
        seg.hcnwino.k2 = 1e-4 # 1e-5 # 
        seg.hcnwino.ghbar = gbar_h
      # Na
      sec.gbar_nax = gbar_nax
      sec.sh_nax = sh_nax
      # KDR
      sec.gbar_kdr = kdr_gbar
      # K-A
      sec.gbar_kap = kap_gbar
      sec.sh_kap = sh_kap
    soma = self.soma
    soma.gbar_kdmc  = kdmc_gbar * kdmc_gbar_somam
    soma.gbar_nax = nax_gbar * nax_gbar_somam
    soma.gbar_kdr = kdr_gbar * kdr_gbar_somam
    soma.gbar_kap = kap_gbar * kap_gbar_somam
    soma.gkbar_ikc = ikc_gkbar
    soma.gcalbar_cal = cal_gcalbar * cal_gbar_somam
    soma.gcanbar_can = can_gcanbar * can_gbar_somam
    soma.gcatbar_cat = cat_gcatbar * cat_gbar_somam
    h.distance(0,0.5,sec=self.soma) # middle of soma is origin for distance
    for sec in self.apic:
      sec.gcalbar_cal = cal_gcalbar
      sec.gcanbar_can = can_gcanbar
      sec.gcatbar_cat = cat_gcatbar
      sec.gkbar_ikc = ikc_gkbar * ikc_gbar_dendm
      sec.gbar_cagk = cagk_gbar
      for seg in sec:
        d = h.distance(seg.x,sec=sec)
        seg.hcnwino.ghbar = gbar_h * exp(d/h_lambda)
        seg.gmax_km = km_gmax * exp(d/h_lambda)
        seg.gbar_kap = soma.gbar_kap * exp(d/h_lambda)
        seg.gbar_kdr = soma.gbar_kdr * exp(d/h_lambda)
    self.apic[1].gcalbar_cal = cal_gcalbar * calginc # middle apical dend gets more iL
    self.apic[2].cm = 2.0

    Bdend = self.Bdend
    Bdend.gcalbar_cal = cal_gcalbar * cal_gbar_bdendm
    Bdend.gcanbar_can = can_gcanbar * can_gbar_bdendm
    Bdend.gcatbar_cat = cat_gcatbar * cat_gbar_bdendm
    Bdend.gkbar_ikc = ikc_gkbar * ikc_gbar_dendm
    Bdend.gbar_cagk = cagk_gbar
    Bdend.gbar_kap = soma.gbar_kap; Bdend.gbar_kdr = soma.gbar_kdr
    Bdend.gmax_km = km_gmax

  def set_conductances (self): # insert the conductances
    for sec in self.all_sec:
      sec.insert('k_ion')
      sec.insert('na_ion')
      sec.insert('ca_ion')
      sec.insert('pas')         # passive   
      sec.insert('hcnwino') # H channel in Ih.mod
      sec.insert('nax')      # Na current
      sec.insert('kdr')      # K delayed rectifier current
      sec.insert('kap')      # K-A current
      # calcium-related channels
      sec.insert('cal') # cal_mig.mod
      sec.insert('can') # can_mig.mod
      sec.insert('cat') # cat_mig.mod
      sec.insert('ikc') # IC.mod - ca and v dependent k channel - BK
    soma = self.soma; self.soma.insert('kdmc')  # K-D current in soma only
    for sec in self.apic:
      sec.insert('km') # km.mod
      sec.insert('cagk') # cagk.mod - SK
    self.Bdend.insert('km') # km.mod
    self.Bdend.insert('cagk') # cagk.mod - SK
		
  def set_synapses(self):
    erevgaba = -80
    self.somaGABAf = Synapse(sect=self.soma,loc=0.5,tau1=0.07,tau2=9.1,e=erevgaba)
    self.somaAMPA = Synapse(sect=self.soma,loc=0.5,tau1=0.05,tau2=5.3,e=0)
    bdsyloc = 0.5 
    self.BdendAMPA = Synapse(sect=self.Bdend,loc=bdsyloc,tau1=0.05, tau2=5.3,e=0)    
    self.BdendNMDA = SynapseNMDA(sect=self.Bdend,loc=bdsyloc,tau1NMDA=tau1NMDAEE,tau2NMDA=tau2NMDAEE,r=1,e=0)
    self.Adend1GABAs = Synapse(sect=self.Adend1,loc=0.5,tau1=0.2,tau2=20,e=erevgaba)
    self.Adend2GABAs = Synapse(sect=self.Adend2,loc=0.5,tau1=0.2,tau2=20,e=erevgaba)
    self.Adend3GABAs = Synapse(sect=self.Adend3,loc=0.5,tau1=0.2,tau2=20,e=erevgaba)
    self.Adend3GABAf = Synapse(sect=self.Adend3,loc=0.5,tau1=0.07,tau2=9.1,e=erevgaba)
    self.Adend3AMPA = Synapse(sect=self.Adend3,loc=0.5,tau1=0.05,tau2=5.3,e=0)
    self.Adend3NMDA = SynapseNMDA(sect=self.Adend3,loc=0.5,tau1NMDA=tau1NMDAEE,tau2NMDA=tau2NMDAEE,r=1,e=0)
    self.Adend2AMPA = Synapse(sect=self.Adend2,loc=0.5,tau1=0.05,tau2=5.3,e=0)
    self.Adend2NMDA = SynapseNMDA(sect=self.Adend2,loc=0.5,tau1NMDA=tau1NMDAEE,tau2NMDA=tau2NMDAEE,r=1,e=0)
    self.Adend1AMPA = Synapse(sect=self.Adend1,loc=0.5,tau1=0.05,tau2=5.3,e=0)
    self.Adend1NMDA = SynapseNMDA(sect=self.Adend1,loc=0.5,tau1NMDA=tau1NMDAEE,tau2NMDA=tau2NMDAEE,r=1,e=0)
    self.Adend3mGLUR = SynapsemGLUR(sect=self.Adend3,loc=0.5)
    self.Adend3GABAB = SynapseGABAB(sect=self.Adend3,loc=0.5)
    self.Adend2mGLUR = SynapsemGLUR(sect=self.Adend2,loc=0.5)
    self.Adend2GABAB = SynapseGABAB(sect=self.Adend2,loc=0.5)
    self.Adend1mGLUR = SynapsemGLUR(sect=self.Adend1,loc=0.5)
    self.Adend1GABAB = SynapseGABAB(sect=self.Adend1,loc=0.5)

#######################################
#      some utils to avoid the h.     #
vlk = h.vlk
Vector = h.Vector
NQS = h.NQS
gg = h.gg
ge = h.ge
Random = h.Random
List = h.List
Matrix = h.Matrix
nqsdel = h.nqsdel
Graph = h.Graph
vrsz = h.vrsz
allocvecs = h.allocvecs
NetCon = h.NetCon
NetStim = h.NetStim
#######################################

