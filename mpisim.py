from neuron import *
# MPI
pc = h.ParallelContext() # MPI: Initialize the ParallelContext class
nhosts = int(pc.nhost()) # Find number of hosts
if pc.id()==0: print('\nSetting up network...')
import sys
import os
import string
h("strdef simname, allfiles, simfiles, output_file, datestr, uname, osname, comment")
h.simname=simname = "dystdemo"
h.allfiles=allfiles = "pyinit.py geom.py mpisim.py"
h.simfiles=simfiles = "pyinit.py geom.py mpisim.py"
h("runnum=1")
runnum = 1.0
h.datestr=datestr = "16jun21"
h.output_file=output_file = "data/16jun21.14"
h.uname=uname = "x86_64"
h.osname=osname="linux"
h("templates_loaded=0")
templates_loaded=0
h("xwindows=1.0")
xwindows = 1.0
h.xopen("nrnoc.hoc")
h.xopen("init.hoc")
h("proc setMemb () { }") # so e_pas will not get modified
CTYPi = 60.0
STYPi = 20.0
from pyinit import *
from labels import *
delm = numpy.zeros( (CTYPi, CTYPi) )
deld = numpy.zeros( (CTYPi, CTYPi) )
pmat = numpy.zeros( (CTYPi, CTYPi) )
synloc = numpy.zeros( (CTYPi, CTYPi) )
from geom import *
from vector import *
from nqs import *
import random
from pylab import *
from datetime import datetime

#########################################################################
# global params
verbose = dconf['verbose']
ISEED = dconf['iseed']
WSEED = dconf['wseed']
PSEED = dconf['pseed']
scale = dconf['scale']
gGID = 0 # global ID for cells
pmatscale = 1.0 # 1.0 / scale
spiketh = -15 # spike threshold, 10 mV is NetCon default, lower it for all cells

simstr = dconf['simstr']
saveout = dconf['saveout']
recdt = dconf['recdt']
recvdt = dconf['recvdt']
saveconns = dconf['saveconns']
indir = dconf['indir']
outdir = dconf['outdir']

# make dir, catch exceptions
def safemkdir (dn):
  try:
    os.mkdir(dn)
    return True
  except OSError:
    if not os.path.exists(dn):
      print 'could not create', dn
      return False
    else:
      return True

# backup the config file
def backupcfg (simstr):
  safemkdir('backupcfg')
  fout = 'backupcfg/' + simstr + '.cfg'
  if os.path.exists(fout):
    print 'removing prior cfg file' , fout
    os.system('rm ' + fout)  
  os.system('cp ' + fcfg + ' ' + fout) # fcfg created in geom.py via conf.py

if pc.id()==0: 
  backupcfg(simstr) # backup the config file
  print "config file is " , fcfg

h.tstop = tstop = dconf['tstop']
tstart = 0.0
h.dt = dconf['dt']
h.steps_per_ms = 1/h.dt
h.v_init = -65
h.celsius = 37
h.fracca_MyExp2SynNMDABB = dconf['nmfracca'] # fraction of NMDA current that is from calcium
rdmsec = dconf['rdmsec']

EEGain = dconf['EEGain']
EIGainFS = dconf['EIGainFS']
EIGainLTS = dconf['EIGainLTS']
IEGain = dconf['IEGain']
IIGain = dconf['IIGain']
IIGainLTSFS =  IIGain 
IIGainFSLTS =  IIGain 
IIGainLTSLTS = IIGain 
IIGainFSFS =   IIGain 
GB2R = dconf['GB2R']; 
NMAMREE = dconf['NMAMREE'] 
NMAMREI = dconf['NMAMREI'] 
mGLURR = dconf['mGLURR'] # ratio of mGLUR weights to AM2 weights
cpernet = []  # cells of a given type for network
wmat = numpy.zeros( (CTYPi, CTYPi, STYPi) ) # internal weights
wmatex = numpy.zeros( (CTYPi, STYPi) ) # external weights
ratex = numpy.zeros( (CTYPi, STYPi) )  # external rates
EXGain = dconf['EXGain']
sgrhzEE = dconf['sgrhzEE'] # external E inputs to E cells; 1000 is default
sgrhzEI = dconf['sgrhzEI'] # external E inputs to I cells
sgrhzIE = dconf['sgrhzIE'] # external I inputs to E cells
sgrhzII = dconf['sgrhzII'] # external I inputs to I cells
sgrhzNME = dconf['sgrhzNME'] # external NM inputs to E cells; 10 is default
sgrhzNMI = dconf['sgrhzNMI'] # external NM inputs to I cells
sgrhzMGLURE = dconf['sgrhzMGLURE'] # external mGLUR inputs to E cells
sgrhzGB2 = dconf['sgrhzGB2'] # external inputs onto E cell GB2 synapses

# params for swire
colside = 120 # 120 microns squared
slambda = 100 # space constant for wiring falloff (used for I -> X)
axonalvelocity = 10000 # axonal velocity in um/ms -- this is 10 mm/s
#########################################################################

# setwmatex - set weights of external inputs to cells
def setwmatex ():
  for ct in xrange(CTYPi):
    for sy in xrange(STYPi):
      ratex[ct][sy]=0
      wmatex[ct][sy]=0
  for ct in xrange(CTYPi):
    if cpernet[ct] <= 0: continue
    if IsLTS(ct): # dendrite-targeting interneurons (LTS cells)
      ratex[ct][AM2]=sgrhzEI
      ratex[ct][NM2] = sgrhzNMI
      ratex[ct][GA]=sgrhzII
      ratex[ct][GA2]=sgrhzII
      wmatex[ct][AM2] = 0.02e-3 
      wmatex[ct][NM2] = 0.02e-3
      wmatex[ct][GA]=  0
      wmatex[ct][GA2]= 0.2e-3 # * 0
    elif ice(ct): # soma-targeting interneurons (basket/FS cells)
      ratex[ct][AM2]=sgrhzEI
      ratex[ct][NM2] = sgrhzNMI
      ratex[ct][GA]=sgrhzII
      ratex[ct][GA2]=sgrhzII
      wmatex[ct][AM2] = 0.02e-3 * 5.0
      wmatex[ct][NM2] = 0.02e-3 * 5.0
      wmatex[ct][GA]= 0
      wmatex[ct][GA2]= 0.2e-3 
    else: # E cells
      ratex[ct][MG]=sgrhzMGLURE 
      ratex[ct][AM2]=sgrhzEE
      ratex[ct][NM2]=sgrhzNME
      ratex[ct][GA]=sgrhzIE
      ratex[ct][GA2]=sgrhzIE
      ratex[ct][GB2]=sgrhzGB2
      wmatex[ct][MG] = 1.5 
      wmatex[ct][AM2] = 0.02e-3
      wmatex[ct][NM2] = 0.02e-3
      wmatex[ct][GA] = 0.2e-3
      wmatex[ct][GA2] = 0.2e-3 
      wmatex[ct][GB2] = 5e-3
    for sy in xrange(STYPi): wmatex[ct][sy] *= EXGain # apply gain control

# set number of cells of a type in the network at scale==1
def setcpernet ():
  global cpernet
  cpernet = []
  for i in xrange(CTYPi): cpernet.append(0)
  cpernet[E2]  = 300
  cpernet[I2]  =  37
  cpernet[I2L] =  37
  cpernet[E4] =  173
  cpernet[E5R] = 150
  cpernet[E5B] = 196
  cpernet[E5P] = 196
  cpernet[I5]  =  89
  cpernet[I5L] =  89
  cpernet[E6] =  179
  cpernet[E6C] = 179
  cpernet[I6] =   45
  cpernet[I6L] =  45

# synapse locations DEND SOMA AXON
def setsynloc ():
  for ty1 in xrange(CTYPi):
    for ty2 in xrange(CTYPi):
      if ice(ty1):
        if IsLTS(ty1):
          synloc[ty1][ty2]=DEND # distal [GA2] - from LTS
        else:
          synloc[ty1][ty2]=SOMA # proximal [GA] - from FS
      else:
        synloc[ty1][ty2]=DEND # E always distal. use AM2,NM2

# setdelmats -- setup delm,deld
def setdelmats ():
  for ty1 in xrange(CTYPi):
    for ty2 in xrange(CTYPi):
      if synloc[ty1][ty2]==DEND and ice(ty2):
        # longer delays at dendrites of interneurons since they are currently single compartment
        delm[ty1][ty2]=2.0 
        deld[ty1][ty2]=0.2 
      else:
        delm[ty1][ty2]=2.0
        deld[ty1][ty2]=0.2

# weight params
def setwmat ():
  for ty1 in xrange(CTYPi):
    for ty2 in xrange(CTYPi):
      for sy in xrange(STYPi): wmat[ty1][ty2][sy]=0
  # E2 -> I weight
  wmat[E2][I2][AM2] = wmat[E2][I2L][AM2] = 0.78
  wmat[E2][I5][AM2] = 0.11
  wmat[E2][I5L][AM2] = 1.01
  # E4 -> I weight
  wmat[E4][I2][AM2] = wmat[E4][I2L][AM2] = 0.3625
  wmat[E4][I5][AM2] = 1.0775
  wmat[E4][I5L][AM2] = 0.1225
  wmat[E4][I6][AM2] = wmat[E4][I6L][AM2] = 0.4375
  # E5R (IT E5a) -> I weight
  wmat[E5R][I2][AM2] = wmat[E5R][I2L][AM2] = 0.3625
  wmat[E5R][I5][AM2] = 1.0775
  wmat[E5R][I5L][AM2] = 0.1225
  wmat[E5R][I6][AM2] = wmat[E5R][I6L][AM2] = 0.4375
  # E5B (IT E5b) -> I weight
  wmat[E5B][I2][AM2] = wmat[E5B][I2L][AM2] = 0.3625
  wmat[E5B][I5][AM2] = 1.0775
  wmat[E5B][I5L][AM2] = 0.1225
  wmat[E5B][I6][AM2] = wmat[E5B][I6L][AM2] = 0.4375
  # E5P (PT E5b) -> I weight
  wmat[E5P][I2][AM2] = wmat[E5P][I2L][AM2] = 0.3625
  wmat[E5P][I5][AM2] = 1.0775
  wmat[E5P][I5L][AM2] = 0.1225
  wmat[E5P][I6][AM2] = wmat[E5P][I6L][AM2] = 0.4375
  # E6 (IT) -> I weight
  wmat[E6][I5][AM2] = wmat[E6][I5L][AM2] = 0.24786
  wmat[E6][I6][AM2] = wmat[E6][I6L][AM2] = 0.53
  # E6C (CT) -> I weight
  wmat[E6C][I5][AM2] = wmat[E6C][I5L][AM2] = 0.24786
  wmat[E6C][I6][AM2] = wmat[E6C][I6L][AM2] = 0.53
  epops = [E2, E4, E5R, E5B, E5P, E6, E6C]
  prty = E2; # E2 -> E weight
  for poty,wght in zip(epops,[0.6416, 0.3700, 0.5049, 0.4446, 0.4446, 0.0000, 0.0000]): wmat[prty][poty][AM2] = wght  
  prty = E4  # E4 -> E weight
  for poty,wght in zip(epops,[0.7368, 0.6416, 0.6383, 0.8981, 0.8981, 1.9082, 1.9082]): wmat[prty][poty][AM2] = wght
  prty = E5R # E5R -> E weight
  for poty,wght in zip(epops,[0.5243, 0.4180, 0.5716, 0.8320, 0.8320, 0.3220, 0.3220]): wmat[prty][poty][AM2] = wght
  prty = E5B # E5B -> E weight
  for poty,wght in zip(epops,[0.2342, 0.1700, 0.3194, 0.6855, 0.6855, 0.4900, 0.4900]): wmat[prty][poty][AM2] = wght
  # E5P -> E weight
  wmat[E5P][E5P][AM2] = 0.6855
  prty = E6 # E6 -> E weight
  for poty,wght in zip(epops,[0.0000, 0.0000, 0.1365, 0.3092, 0.3092, 0.5300, 0.5300]): wmat[prty][poty][AM2] = wght
  prty = E6C # E6C -> E weight
  for poty,wght in zip(epops,[0.0000, 0.0000, 0.1365, 0.3092, 0.3092, 0.5300, 0.5300]): wmat[prty][poty][AM2] = wght
  # I -> X
  for prty,sy in zip([I2,I2L],[GA,GA2]):
    for poty in [E2,I2,I2L]: 
      wmat[prty][poty][sy] = 1.5
      if IsLTS(prty) and not ice(poty): wmat[prty][poty][GB2] = 1.5 * GB2R
  for prty,sy in zip([I5,I5L],[GA,GA2]):
    for poty in [E4,E5R,E5B,E5P,I5,I5L]: 
      wmat[prty][poty][sy] = 1.5
      if IsLTS(prty) and not ice(poty): wmat[prty][poty][GB2] = 1.5 * GB2R
  for prty,sy in zip([I6,I6L],[GA,GA2]):
    for poty in [E6,E6C,I6,I6L]: 
      wmat[prty][poty][sy] = 1.5
      if IsLTS(prty) and not ice(poty): wmat[prty][poty][GB2] = 1.5 * GB2R
  # gain control
  for ty1 in xrange(CTYPi):
    for ty2 in xrange(CTYPi):
      for sy in xrange(STYPi):
        if wmat[ty1][ty2][sy] > 0:
          if ice(ty1): # I -> X
            if ice(ty2):
              if IsLTS(ty1): # LTS -> I
                if IsLTS(ty2): # LTS -> LTS
                  gn = IIGainLTSLTS
                else: # LTS -> FS
                  gn = IIGainLTSFS
              else: # FS -> I
                if IsLTS(ty2): # FS -> LTS
                  gn = IIGainFSLTS
                else: # FS -> FS
                  gn = IIGainFSFS
            else: # I -> E
              gn = IEGain
          else: # E -> X
            if ice(ty2): # E -> I
              if IsLTS(ty2): # E -> LTS
                gn = EIGainLTS
              else: # E -> FS
                gn = EIGainFS
            else: # E -> E
              gn = EEGain
              if sy==AM2: 
                wmat[ty1][ty2][MG] = wmat[ty1][ty2][AM2] * mGLURR
                if verbose: print 'AM2:',wmat[ty1][ty2][AM2],'mGLURR:',wmat[ty1][ty2][MG]
            if sy==AM2:
              if ice(ty2): # E -> I
                wmat[ty1][ty2][NM2] = wmat[ty1][ty2][AM2] * NMAMREI
              else: # E -> E
                wmat[ty1][ty2][NM2] = wmat[ty1][ty2][AM2] * NMAMREE
          wmat[ty1][ty2][sy] *= gn 

def setpmat ():
  for ii in xrange(CTYPi):
    for jj in xrange(CTYPi): pmat[ii][jj]=0
  # E2 -> I wiring
  pmat[E2][I2] = pmat[E2][I2L] = 0.1871
  pmat[E2][I5] = 0.02
  pmat[E2][I5L] = 0.217
  # E4 -> I wiring
  pmat[E4][I2] = pmat[E4][I2L] = 0.0222
  pmat[E4][I5] = 0.1906
  pmat[E4][I5L] = 0.0349
  pmat[E4][I6] = pmat[E4][I6L] = 0.0155
  # E5R (IT E5a) -> I wiring
  pmat[E5R][I2] = pmat[E5R][I2L] = 0.0222
  pmat[E5R][I5] = 0.1906
  pmat[E5R][I5L] = 0.0349
  pmat[E5R][I6] = pmat[E5R][I6L] = 0.0155
  # E5B (IT E5b) -> I wiring
  pmat[E5B][I2] = pmat[E5B][I2L] = 0.0222
  pmat[E5B][I5] = 0.1906
  pmat[E5B][I5L] = 0.0349
  pmat[E5B][I6] = pmat[E5B][I6L] = 0.0155
  # E5P (PT E5b) -> I wiring
  pmat[E5P][I2] = pmat[E5P][I2L] = 0.0222
  pmat[E5P][I5] = 0.1906
  pmat[E5P][I5L] = 0.0349
  pmat[E5P][I6] = pmat[E5P][I6L] = 0.0155
  # E6 (IT) -> I wiring
  pmat[E6][I5] = pmat[E6][I5L] = 0.0249
  pmat[E6][I6] = pmat[E6][I6L] = 0.0234
  # E6C (CT) -> I wiring
  pmat[E6C][I5] = pmat[E6C][I5L] = 0.0249
  pmat[E6C][I6] = pmat[E6C][I6L] = 0.0234
  epops = [E2, E4, E5R, E5B, E5P, E6, E6C]  
  prty = E2;   # E2 -> E wiring
  for poty,prob in zip(epops,[0.1503, 0.1100, 0.0509, 0.0124, 0.0690, 0,0]): pmat[prty][poty] = prob    
  prty = E4    # E4 -> E wiring
  for poty,prob in zip(epops,[0.0523, 0.1503, 0.0413, 0.0143, 0.0120, 0.0030, 0.0030]): pmat[prty][poty] = prob  
  prty = E5R   # E5R -> E wiring
  for poty,prob in zip(epops,[0.0355, 0.0275, 0.1810, 0.0142, 0.0205, 0.0136, 0.0136]): pmat[prty][poty] = prob
  prty = E5B   # E5B -> E wiring
  for poty,prob in zip(epops,[0.0167, 0.0293, 0.0536, 0.1810, 0.0448, 0.0160, 0.0160]): pmat[prty][poty] = prob  
  pmat[E5P][E5P] = 0.1810 # E5P -> E5P  
  prty = E6 # E6 -> E wiring
  for poty,prob in zip(epops,[0, 0, 0.0334, 0.0271, 0.0277, 0.0282, 0.0234]): pmat[prty][poty] = prob  
  prty = E6C # E6C -> E wiring
  for poty,prob in zip(epops,[0, 0, 0.0334, 0.0271, 0.0277, 0.0234, 0.0282]): pmat[prty][poty] = prob
  # I -> X
  for prty in [I2,I2L]:
    for poty in [E2,I2,I2L]: pmat[prty][poty] = 1.0
  for prty in [I5,I5L]:
    for poty in [E4,E5R,E5B,E5P,I5,I5L]: pmat[prty][poty] = 1.0
  for prty in [I6,I6L]:
    for poty in [E6,E6C,I6,I6L]: pmat[prty][poty] = 1.0  
  for ii in xrange(CTYPi):
    for jj in xrange(CTYPi): pmat[ii][jj]*=pmatscale

numc = [0 for i in xrange(CTYPi)]; # number of cells of a type
ix = [0 for i in xrange(CTYPi)]; #starting index of a cell type (into self.ce list)
ixe = [0 for i in xrange(CTYPi)]; #ending index of a cell type
allcells,ecells,icells = 0,0,0
div = zeros( (CTYPi, CTYPi) )
conv = zeros( (CTYPi, CTYPi) )
syty1 = zeros( (CTYPi, CTYPi) ) # stores synapse codes (from labels.py)
syty2 = zeros( (CTYPi, CTYPi) ) # stores synapse code (from labels.py)
syty3 = zeros( (CTYPi, CTYPi) ) # stores synapse code (from labels.py)
sytys1 = {} # dictionary of synapse names
sytys2 = {} # dictionary of synapse names
sytys3 = {} # dictionary of synapse names
SOMA = 0; BDEND = 1; ADEND1 = 2; ADEND2 = 3; ADEND3 = 4;
dsecnames = ['soma','Bdend','Adend1','Adend2','Adend3']

def setdivmat ():
  import math
  for ty1 in xrange(CTYPi):
    for ty2 in xrange(CTYPi):
      if pmat[ty1][ty2] > 0.0: 
        div[ty1][ty2] =  math.ceil(pmat[ty1][ty2]*numc[ty2])
        conv[ty1][ty2] = int(0.5 + pmat[ty1][ty2]*numc[ty1])

# setup cell-type-to-cell-type synapse-type information
def setsyty ():
  for ty1 in xrange(CTYPi): # go thru presynaptic types
    for ty2 in xrange(CTYPi): # go thru postsynaptic types
      syty1[ty1][ty2] = syty2[ty1][ty2] = syty3[ty1][ty2] = -1 # initialize to invalid
      if numc[ty1] <= 0 or numc[ty2] <= 0: continue
      if ice(ty1): # is presynaptic type inhibitory?
        if IsLTS(ty1): # LTS -> X
          syty1[ty1][ty2] = GA2 # code for dendritic gabaa synapse
          if ice(ty2): # LTS -> Io
            sytys1[(ty1,ty2)] = "GABAss"
          else: # LTS -> E
            syty2[ty1][ty2] = GB2 # code for denritic gabab synapse
            sytys1[(ty1,ty2)] = "GABAs"
            sytys2[(ty1,ty2)] = "GABAB"
        else: # BAS -> X
          syty1[ty1][ty2] = GA # code for somatic gabaa synapse
          sytys1[(ty1,ty2)] = "GABAf"
      else: # E -> X
        syty1[ty1][ty2] = AM2 # code for dendritic ampa synapse
        syty2[ty1][ty2] = NM2 # code for dendritic nmda synapse
        if ice(ty2): # E -> I
          sytys1[(ty1,ty2)] = "AMPA"
          sytys2[(ty1,ty2)] = "NMDA"
        else: # E -> E
          sytys1[(ty1,ty2)] = "AMPA"
          sytys2[(ty1,ty2)] = "NMDA"
          sytys3[(ty1,ty2)] = "mGLUR"
          syty3[ty1][ty2] = MG # use MG -- for mGluR

lctyID,lctyClass = [],[]

# setup some convenient data structures
def setix (scale):
  import math
  global allcells,ecells,icells
  for i in xrange(CTYPi):
    numc[i] = int(math.ceil(cpernet[i]*scale))
    if numc[i] > 0:
      ty = PyrAdr
      if ice(i):
        if IsLTS(i): ty = LTS
        else: ty = FS
      for j in xrange(numc[i]):
        lctyClass.append(ty)
        lctyID.append(i)
      allcells += numc[i]
      if ice(i): icells += numc[i]
      else: ecells += numc[i]
  sidx = 0
  for i in xrange(CTYPi):
    if numc[i] > 0:
      ix[i] = sidx
      ixe[i] = ix[i] + numc[i] - 1
      sidx = ixe[i] + 1
  setdivmat()
  setsyty()

# setcellpos([pseed,network diameter in microns])
def setcellpos (pseed=4321,cside=colside):
  rdm=Random(); rdm.ACG(pseed)
  cellsnq = NQS("id","ty","ice","xloc","yloc","zloc")
  cellsnq.clear(allcells) # alloc space
  lX,lY,lZ=[],[],[]
  ldy = {}
  ldy[E2] =  (160, 420)
  ldy[E4] = (420, 570)
  ldy[E5R] = (570, 700)
  ldy[E5B] = (700, 1040)
  ldy[E5P] = (700, 1040)
  ldy[E6] =  (1040, 1350)
  ldy[E6C] = (1040, 1350)
  ldy[I2] = (160, 420)
  ldy[I2L] = (160, 420)
  ldy[I5] = (420, 1040)
  ldy[I5L] = (420, 1040)
  ldy[I6] = (1040, 1350)
  ldy[I6L] = (1040, 1350)
  for i in xrange(allcells):    
    ctyp = lctyID[i]
    [x,y,z] = [rdm.uniform(0,cside), rdm.uniform(ldy[ctyp][0],ldy[ctyp][1]), rdm.uniform(0,cside)]
    cellsnq.append(i,ctyp,ice(ctyp),x,y,z); lX.append(x); lY.append(y); lZ.append(z);
  return cellsnq,lX,lY,lZ

setcpernet() # setup number of cells per network
setwmatex() # setup matrices of external inputs
setsynloc() # setup synapse location matrices
setdelmats() # setup delay matrices
setwmat() # setup weight matrix
setpmat() # setup connectivity matrix
setix(scale)
cellsnq,lX,lY,lZ=setcellpos()

ce = [] # cells on the host
gidvec = [] # gids of cells on the host
lncrec,ltimevec,lidvec=[],[],[] # spike recorders
dlids = {} # map from gid back to ce index

# create the cells
pcID = int(pc.id()); 
maxcells=0
cperhost = int(allcells/nhosts)
maxcells = cperhost
extra = allcells - cperhost*nhosts
if extra > 0: # check if any remainder cells
  if pcID < extra: # first hosts get extra cell
    maxcells += 1 # assign an extra cell if any remainders
    gid = pcID * (cperhost + 1)
  else: # rest of hosts do not
    gid = extra*(cperhost+1) + (pcID-extra) * cperhost
else: # even division? all hosts get equal cells
  gid = pcID * cperhost
for i in xrange(maxcells):
  ct = lctyID[gid]
  cell = lctyClass[gid](0+i*50,0,0,gid,ct)
  cell.x,cell.y,cell.z = lX[gid],lY[gid],lZ[gid]
  dlids[gid] = len(ce) # map from gid back to ce index
  ce.append(cell)
  gidvec.append(gid)
  pc.set_gid2node(gid,pcID)
  timevec,idvec = h.Vector(),h.Vector()
  ncrec = h.NetCon(ce[-1].soma(0.5)._ref_v, None, sec=ce[-1].soma)
  ncrec.record(timevec,idvec,gid)
  ncrec.threshold = spiketh # 10 mV is default, lower it for FS cells
  ltimevec.append(timevec); lidvec.append(idvec); lncrec.append(ncrec)
  pc.cell(gid,lncrec[-1],1) # 1 as 3rd arg means this cell can be source for events
  gid += 1
  
print('  Number of cells on node %i: %i' % (pcID,len(ce)))
pc.barrier()

# wire the network using a NQS table
nccl = []
def wirenq (cnq):
  global nccl
  nccl = [] # NetCon list for connections between cells 
  cnq.tog("DB")
  vid1,vid2,vwt1,vwt2,vdel,vsec=cnq.getcol("id1"),cnq.getcol("id2"),cnq.getcol("wt1"),cnq.getcol("wt2"),cnq.getcol("del"),cnq.getcol("sec")
  vwt3 = cnq.getcol("wt3")
  for i in xrange(int(cnq.v[0].size())):
    prid = int(vid1[i])
    poid = int(vid2[i])
    if not pc.gid_exists(poid): continue # only make the connection on a node that has the target
    ty1 = lctyID[prid] 
    ty2 = lctyID[poid] 
    sname = dsecnames[int(vsec[i])] # which section is the synapse on?
    syn = sname + sytys1[(ty1,ty2)]
    wt1 = vwt1[i]
    delay = vdel[i]
    targ = ce[dlids[poid]]
    nc1 = pc.gid_connect(prid, targ.__dict__[syn].syn)
    nc1.delay = delay; nc1.weight[0] = wt1; nc1.threshold = spiketh; nccl.append(nc1)
    wt2 = vwt2[i]
    if wt2 > 0: # two synapses? (i.e., AMPA and NMDA)
      syn = sname + sytys2[(ty1,ty2)]
      if syn in targ.__dict__: # since GABAB not in Bdend
        nc2 = pc.gid_connect(prid, targ.__dict__[syn].syn)
        nc2.delay = delay; nc2.weight[0] = wt2; nc2.threshold = spiketh; nccl.append(nc2)
    wt3 = vwt3[i]
    if wt3 > 0: # three synapses? (i.e., AMPA and NMDA and mGLUR)
      if verbose: print 'mGLUR synapse wt3 > 0:',wt3
      syn = sname + sytys3[(ty1,ty2)]
      if syn in targ.__dict__: # make sure target has this synapse (needed since Bdend does not have mGLUR)
        nc3 = pc.gid_connect(prid, targ.__dict__[syn].syn)
        nc3.delay = delay; nc3.weight[0] = wt3; nc3.threshold = spiketh; nccl.append(nc3)

#
def picksec (prty, poty, rdm):
  if ice(poty): return SOMA
  if ice(prty): # I -> E
    if IsLTS(prty): # LTS -> E
      if rdmsec: return rdm.discunif(ADEND1,ADEND3)
      else: return ADEND3
    else:
      return SOMA
  else: # E -> E
    if rdmsec: return rdm.discunif(BDEND,ADEND3)
    else: return ADEND3

# swire - spatial wiring: wires the network using pmat and cell positions
#                    (wiring probability affected by distance btwn cells)
#  slambda (global) specifies length-constant for spatially-dependent fall-off in wiring probability
#  slambda is only used for I->X wiring; for E->X wiring uses fixed probability
def swire (wseed):
  global slambda
  from math import sqrt,exp
  [vidx,vdel,vtmp,vwt1,vwt2,vwt3,vprob] = [Vector() for x in xrange(7)]
  z = 0
  if slambda <= 0:
    print "swire WARN: invalid slambda=", slambda, "setting slambda to ", colside/3
    slambda=colside/3
  slambdasq = slambda**2 # using squared distance
  h.vrsz(1e4,vidx,vdel,vtmp)
  rdm=Random(); rdm.ACG(wseed) #initialize random # generator
  rdm.uniform(0,1)
  vprob.resize(allcells**2); vprob.setrand(rdm)
  pdx=0 # index into vprob
  connsnq=NQS("id1","id2","del","wt1","wt2","wt3","sec")
  connsnq.clear(1e3*allcells)
  for prid in xrange(allcells): 
    vrsz(0,vidx,vdel,vwt1,vwt2,vwt3)
    prty=lctyID[prid]
    ic1=ice(prty)
    for poty in xrange(0,CTYPi):
      if numc[poty] > 0 and pmat[prty][poty]>0:
        pbase = pmat[prty][poty]
        for poid in xrange(ix[poty],ixe[poty]+1): # go thru postsynaptic cells
          if prid==poid: continue # no self-connects
          ic2=ice(lctyID[poid])
          dx = lX[prid] - lX[poid]
          dy = lY[prid] - lY[poid]
          dz = lZ[prid] - lZ[poid]
          ds = sqrt(dx**2 + dy**2 + dz**2) # Connectivity fall-off depends on 3D distance
          if ic1: prob = exp(-ds/slambda) # probability of connect falls off only for I->X
          else: prob = pbase
          if prob >= vprob[pdx]: # rdm.uniform(0,1)
            mindelay = delm[prty][poty]-deld[prty][poty]
            maxdelay = delm[prty][poty]+deld[prty][poty]
            delay=rdm.uniform(mindelay,maxdelay) # synaptic delay
            delay += ds/axonalvelocity # add axonal delay 
            vidx.append(poid); vdel.append(delay)
            if syty1[prty][poty]>=0: vwt1.append(wmat[prty][poty][int(syty1[prty][poty])])
            else: vwt1.append(0)
            if syty2[prty][poty]>=0: vwt2.append(wmat[prty][poty][int(syty2[prty][poty])])
            else: vwt2.append(0)
            if syty3[prty][poty]>=0: vwt3.append(wmat[prty][poty][int(syty3[prty][poty])])
            else: vwt3.append(0)
          pdx += 1
    for ii in xrange(int(vidx.size())): connsnq.append(prid,vidx[ii],vdel[ii],vwt1[ii],vwt2[ii],vwt3[ii],picksec(prty , lctyID[int(vidx[ii])], rdm))
  wirenq(connsnq) # do the actual wiring based on self.connsnq
  return connsnq

connsnq=swire(WSEED)

pc.barrier() # wait for wiring to get completed

# setup rxd for E cells
# get list of all Sections associated with an excitatory cell
def getesec ():
  esec = []
  for cell in ce:
    if ice(cell.ty): continue
    for s in cell.all_sec: esec.append(s)
  return esec
  
def pcidpr (s): 
  global pcID
  print 'host',pcID,':',s

### RXD ###
[cyt,er,cyt_er_membrane,ca,caextrude,serca,leak,CB,caCB,buffering]=[None for i in xrange(10)]
rxdsec=getesec() # Section list for use with rxd
pc.barrier()
if len(rxdsec) > 0: # only use rxd if there are viable Sections
  from neuron import rxd
  rxd.options.use_reaction_contribution_to_jacobian = False # faster (checked a few days before 10/16/13)
  fc, fe = 0.83, 0.17 # cytoplasmic, er volume fractions
  cyt = rxd.Region(rxdsec, nrn_region='i', geometry=rxd.FractionalVolume(fc, surface_fraction=1))
  er  = rxd.Region(rxdsec, geometry=rxd.FractionalVolume(fe))
  cyt_er_membrane = rxd.Region(rxdsec, geometry=rxd.ScalableBorder(1))
  caDiff = 0.233
  ca = rxd.Species([cyt, er], d=caDiff, name='ca', charge=2, initial=dconf['cacytinit'])
  caexinit = dconf['caexinit']
  caextrude = rxd.Rate(ca, (caexinit-ca[cyt])/taurcada, regions=cyt, membrane_flux=False)
  ip3 = rxd.Species(cyt, d=0.283, name='ip3', initial=0.0)
  # action of IP3 receptor
  Kip3=0.13; Kact=0.4
  minf = ip3[cyt] * 1000. * ca[cyt] / (ip3[cyt] + Kip3) / (1000. * ca[cyt] + Kact)
  ip3r_gate_state = rxd.State(cyt_er_membrane, initial=0.8)
  h_gate = ip3r_gate_state[cyt_er_membrane]
  k = dconf['gip3'] * (minf * h_gate) ** 3 
  ip3r = rxd.MultiCompartmentReaction(ca[er]<>ca[cyt], k, k, membrane=cyt_er_membrane)    
  # IP3 receptor gating
  ip3rg = rxd.Rate(h_gate, (1. / (1 + 1000. * ca[cyt] / (0.4)) - h_gate) / 400.0)
  # IP3 degradation - moves towards baseline level (ip3_init)
  ip3degTau = 1000 # 1000 ms
  ip3deg = rxd.Rate(ip3, (0.0-ip3[cyt])/ip3degTau, regions=cyt, membrane_flux=False)
  ### RYR - based on Sneyd et al, 2003
  # constants
  k_a_pos = 1500000000000.0 # mM^-4/ms
  k_a_neg = 0.0288 # /ms
  k_b_pos = 1500000000.0 # mM^-3/ms
  k_b_neg = 0.3859 # /ms
  k_c_pos = 0.00175 # /ms
  k_c_neg = 0.0001 # /ms
  v1ryr = dconf['v1ryr'] # /ms
  Ka_4 = k_a_neg / k_a_pos # Ka**4
  Kb_3 = k_b_neg / k_b_pos # Kb**3
  Kc = k_c_neg / k_c_pos
  # w_state is fraction of RYR not in C2 state (closed state), ie fraction of RYR that is open
  # w_infinity - equ: 29 
  c3ryr = (ca[cyt]**3)/Kb_3; 
  c4ryr = Ka_4/(ca[cyt]**4);
  w_inf = (1.0 + c4ryr + c3ryr) / (1.0+(1.0/Kc)+ c4ryr + c3ryr)
  w_state = rxd.State(cyt_er_membrane, initial=0.9999) # 0)#0.9999) # check if initial == 0???? or can put it as w_inf
  # equ:  8 (which is the same as equ 22)
  w_rate = rxd.Rate(w_state, 1.0 - w_state[cyt_er_membrane] / w_inf)
  # P_ryr - gating variable - equ: 7 (which is the same as 27)
  # (open probability)
  ryr_gate = w_state[cyt_er_membrane] * (1.0 + c3ryr) / (1.0 + c4ryr + c3ryr)
  # the following is extracted from equ 9 and 15
  k_ryr = v1ryr*ryr_gate
  ryr = rxd.MultiCompartmentReaction(ca[er]<>ca[cyt], k_ryr, k_ryr, membrane=cyt_er_membrane)
  def setmGLURflux (): # mGLUR synapses generate ip3 that is fed into rxd ip3
    for c in ce:
      if ice(c.ty): continue
      for syn,seg in zip([c.Adend3mGLUR.syn,c.Adend2mGLUR.syn,c.Adend1mGLUR.syn],[c.Adend3(0.5), c.Adend2(0.5), c.Adend1(0.5)]):
        for node in ip3.nodes(seg): 
          node.include_flux(syn._ref_rip3)
  def setrecip3 ():
    for c in ce:
      if ice(c.ty): continue
      c.soma_ip3cyt = Vector(tstop/h.dt)
      c.soma_ip3cyt.record( ip3[cyt].nodes(c.soma)(0.5)[0]._ref_concentration, recdt )
      c.Adend3_ip3cyt = Vector(tstop/h.dt)
      c.Adend3_ip3cyt.record( ip3[cyt].nodes(c.Adend3)(0.5)[0]._ref_concentration, recdt )    
  # SERCA pump: pumps ca from cyt -> ER
  Kserca = 0.1 # Michaelis constant for SERCA pump
  gserca = dconf['gserca']
  serca = rxd.MultiCompartmentReaction(ca[cyt]>ca[er],gserca*(1e3*ca[cyt])**2/(Kserca**2+(1e3*ca[cyt])**2),membrane=cyt_er_membrane,custom_dynamics=True)
  gleak = dconf['gleak']   # leak channel: bidirectional ca flow btwn cyt <> ER
  leak = rxd.MultiCompartmentReaction(ca[er]<>ca[cyt], gleak, gleak, membrane=cyt_er_membrane)
  def setreccaer (): # setup recording of ca[er] for each pyramidal cell in Adend3,soma center
    for c in ce:
      if ice(c.ty): continue
      c.soma_caer = Vector(tstop/h.dt)
      c.soma_caer.record( ca[er].nodes(c.soma)(0.5)[0]._ref_concentration, recdt )
      c.Adend3_caer = Vector(tstop/h.dt)
      c.Adend3_caer.record( ca[er].nodes(c.Adend3)(0.5)[0]._ref_concentration, recdt )
  CB_init = dconf["CB_init"]
  CB_frate = dconf["CB_frate"]
  CB_brate = dconf["CB_brate"]
  CBDiff = 0.043   # um^2 / msec
  CB = rxd.Species(cyt,d=CBDiff,name='CB',charge=0,initial=CB_init) # CalBindin (Anwar)
  caCB = rxd.Species(cyt,d=CBDiff,name='caCB',charge=0,initial=0.0) # Calcium-CB complex
  kCB = [CB_frate, CB_brate] # forward,backward buffering rates
  buffering = rxd.Reaction(ca+CB <> caCB, kCB[0], kCB[1], regions=cyt)
  def setreccacb (): # setup recording of caCB for each pyramidal cell in Adend3,soma center
    for c in ce:
      if ice(c.ty): continue
      c.soma_caCB = Vector(tstop/h.dt)
      c.soma_caCB.record( caCB.nodes(c.soma)(0.5)[0]._ref_concentration, recdt )
      c.Adend3_caCB = Vector(tstop/h.dt)
      c.Adend3_caCB.record( caCB.nodes(c.Adend3)(0.5)[0]._ref_concentration, recdt )
  setreccaer() # NB: only record from RXD variables after ALL species setup!
  setreccacb() # otherwise, the pointers get messed up.
  setrecip3()
  setmGLURflux()

# setup inputs - first noise inputs
def getsyns ():
  syns = {} # mapping of synapse names, first index is ice, second is synapse code
  syns[ (0,MG) ] = ["Adend3mGLUR","Adend2mGLUR","Adend1mGLUR"]
  syns[ (0,AM2) ] = ["Adend3AMPA","Adend2AMPA","Adend1AMPA","BdendAMPA"]
  syns[ (1,AM2) ] = "somaAMPA"
  syns[ (0,NM2) ] = ["Adend3NMDA","Adend2NMDA","Adend1NMDA","BdendNMDA"]
  syns[ (1,NM2) ] = "somaNMDA"
  syns[ (0,GB2) ] = ["Adend3GABAB","Adend2GABAB","Adend1GABAB"]
  syns[ (0,GA2) ] = ["Adend3GABAs","Adend2GABAs","Adend1GABAs"]
  syns[ (1,GA2) ] = "somaGABAss"
  syns[ (0,GA) ] = "somaGABAf"
  syns[ (1,GA) ] = "somaGABAf"
  return syns

dsstr = ['AMPA', 'NMDA', 'GABAS', 'mGLUR', 'GABAB']

# adds synapses across dendritic fields for the E cells
def addsyns ():
  for cell in ce:
    cell.dsy = {}; cell.vsy = {}
    if ice(cell.ty): continue
    ds = {}; ds[cell.Adend3]='Adend3'; ds[cell.Adend2]='Adend2'; ds[cell.Adend1]='Adend1'; ds[cell.Bdend]='Bdend'
    for sec in [cell.Adend3, cell.Adend2, cell.Adend1, cell.Bdend]:
      llsy = [];
      nloc = sec.nseg
      llvsy = []; # for recording currents
      for i,seg in enumerate(sec):
        if seg.x == 0.0 or seg.x == 1.0: continue # skip endpoints
        lsy = []; loc = seg.x; lvsy = [] #AMPA, NMDA, GABAA_slow, GABAB
        #print 'loc:',loc
        lsy.append(Synapse(sect=sec,loc=loc,tau1=0.05,tau2=5.3,e=0)); lvsy.append(h.Vector())#AMPA
        lsy.append(SynapseNMDA(sect=sec,loc=loc,tau1NMDA=15,tau2NMDA=150,r=1,e=0)) # NMDA
        lvsy.append(h.Vector())
        lsy.append(Synapse(sect=sec,loc=loc,tau1=0.2,tau2=20,e=-80)) # GABAA_slow
        lvsy.append(h.Vector())
        lsy.append(SynapsemGLUR(sect=sec,loc=loc)) # mGLUR
        for node in ip3.nodes(seg): node.include_flux(lsy[-1].syn._ref_rip3 ) # all the sub-segments get flux
        lsy.append(SynapseGABAB(sect=sec,loc=loc)) # GABAB
        lvsy.append(h.Vector())
        llsy.append(lsy); llvsy.append(lvsy)
      cell.dsy[sec] = llsy; cell.vsy[sec] = llvsy
    sec = cell.soma; llsy = []; nloc = sec.nseg; llvsy = []
    for i,seg in enumerate(sec):
      if seg.x == 0.0 or seg.x == 1.0: continue # skip endpoints
      lsy = []; loc = seg.x; lvsy = []
      lsy.append(Synapse(sect=sec,loc=loc,tau1=0.07,tau2=9.1,e=-80)) # GABAA_fast
      lvsy.append(h.Vector())
      lsy.append(Synapse(sect=sec,loc=loc,tau1=0.05,tau2=5.3,e=0) ) # AMPA
      lvsy.append(h.Vector())
      lsy.append(SynapseNMDA(sect=sec,loc=loc,tau1NMDA=15,tau2NMDA=150,r=1,e=0)) # NMDA
      lvsy.append(h.Vector())
      llsy.append(lsy); llvsy.append(lvsy);
    cell.dsy[sec] = llsy; cell.vsy[sec] = llvsy;

addsyns()

#creates NetStims (and associated NetCon,Random) - provide 'noise' inputs
#returns next useable value of sead
def makeNoiseNetStim (cel,nsl,ncl,nrl,nrlsead,syn,w,ISI,time_limit,sead):
  #rd2 = h.Random(); rd2.ACG(sead); rd2.uniform(0,100)
  ns = h.NetStim()
  ns.interval = ISI
  ns.noise = 1			
  ns.number = 2 * time_limit / ISI  # create enough spikes for extra time, in case goes over limit
  if type(syn) == str: nc = h.NetCon(ns,cel.__dict__[syn].syn)
  else: nc = h.NetCon(ns,syn)
  nc.delay = h.dt * 2 # 0
  nc.weight[0] = w
  rds = h.Random()
  rds.negexp(1)            # set random # generator using negexp(1) - avg interval in NetStim
  rds.MCellRan4(sead,sead) # seeds are in order, shouldn't matter			
  ns.noiseFromRandom(rds)  # use random # generator for this NetStim                
  ns.start = tstart # rd2.repick() # start inputs random time btwn 0-1e3 ms to avoid artificial sync
  nsl.append(ns)
  ncl.append(nc)
  nrl.append(rds)
  nrlsead.append(sead)

def makeNoiseNetStims (simdur,rdmseed):
  nsl = [] #NetStim List
  ncl = [] #NetCon List
  nrl = [] #Random List for NetStims
  nrlsead = [] #List of seeds for NetStim randoms
  syns = getsyns() ; 
  for cell in ce: # go through cell types, check weights,rates of inputs
    ct = cell.ty # get cell type code
    if ice(ct): # only has 1 compartment
      for sy in xrange(STYPi):
        if wmatex[ct][sy] <= 0.0 or ratex[ct][sy] <= 0: continue
        syn = syns[(ice(ct),sy)]
        if type(syn) == list:
          for idx,SYN in enumerate(syn):
            makeNoiseNetStim(cell,nsl,ncl,nrl,nrlsead,SYN,wmatex[ct][sy],1e3/ratex[ct][sy],simdur,rdmseed*(cell.ID+1)*(idx+1))
        else:
          makeNoiseNetStim(cell,nsl,ncl,nrl,nrlsead,syn,wmatex[ct][sy],1e3/ratex[ct][sy],simdur,rdmseed*(cell.ID+1))
    else: # E cells - need to distribute noise over all sections
      for sec in [cell.Adend3, cell.Adend2, cell.Adend1]:
        llsy = cell.dsy[sec]
        for lsy in llsy:
          for i,sy in enumerate([AM2,NM2,GA2,MG,GB2]):
            if ratex[ct][sy] > 0. and wmatex[ct][sy] > 0.: 
              makeNoiseNetStim(cell,nsl,ncl,nrl,nrlsead,lsy[i].syn,wmatex[ct][sy],(1e3/ratex[ct][sy]),simdur,rdmseed*(cell.ID+1)*(i+1));
      sec = cell.Bdend; llsy = cell.dsy[sec];
      for lsy in llsy:
        for i,sy in enumerate([AM2,NM2,GA2]):
          if ratex[ct][sy] > 0. and wmatex[ct][sy] > 0.:
            makeNoiseNetStim(cell,nsl,ncl,nrl,nrlsead,lsy[i].syn,wmatex[ct][sy],(1e3/ratex[ct][sy]),simdur,rdmseed*(cell.ID+1)*(i+4)); 
      sec = cell.soma; llsy = cell.dsy[sec];
      for i,sy in enumerate([GA,AM,NM]):
        if ratex[ct][sy] > 0. and wmatex[ct][sy] > 0.:
          for lsy in llsy:
            makeNoiseNetStim(cell,nsl,ncl,nrl,nrlsead,lsy[i].syn,wmatex[ct][sy],(1e3/ratex[ct][sy]),simdur,rdmseed*(cell.ID+1)*(i+7)); rdmseed+=1
  return nsl,ncl,nrl,nrlsead

nsl,ncl,nrl,nrlsead = makeNoiseNetStims(tstart+tstop,ISEED)

pc.barrier() # wait for completion of NetStim creation

#this should be called @ beginning of each sim - done in an FInitializeHandler
def init_NetStims ():
  for i in xrange(len(nrl)):
    rds = nrl[i]
    sead = nrlsead[i]
    rds.MCellRan4(sead,sead)
    rds.negexp(1)			

fihns = h.FInitializeHandler(0, init_NetStims)

# handler for printing out time during simulation run
def fi():
  for i in xrange(int(tstart),int(tstart+tstop),100): h.cvode.event(i, "print " + str(i))

if pc.id() == 0: fih = h.FInitializeHandler(1, fi)

vt=Vector(); vt.record(h._ref_t); # record time

pc.barrier() # wait for NetStims to get setup 

####################################################################################
### simulation run here 
def myrun ():
  pc.set_maxstep(10)
  dastart,daend=None,None
  if pc.id()==0:
    dastart = datetime.now()
    print 'started at:',dastart
  h.stdinit()
  if len(rxdsec)>0: # any sections with rxd?
    ca[er].concentration = dconf['caerinit'] # 100e-6
    ca[cyt].concentration = dconf['cacytinit'] # 100e-6
  pc.psolve(h.t+tstop) # run for tstop
  pc.barrier() # Wait for all hosts to get to this point
  if pc.id()==0:
    daend = datetime.now()
    print 'finished ',tstop,' ms sim at:',daend
    dadiff = daend - dastart;
    print 'runtime:',dadiff, '(',tstop/1e3,' s)'

if dconf['dorun']: myrun()

# concatenate the results so can view/save all at once
lspks,lids=array([]),array([])
for host in xrange(nhosts): # is this loop required? can't just post messages from given host?
  if host == pc.id():
    for i in xrange(len(ltimevec)):
      lspks=concatenate((lspks,ltimevec[i]))
      lids=concatenate((lids,lidvec[i]))    

# save data - output path based on simstr and pcid
def savedata (simstr,pcid):
  safemkdir(outdir)
  fn = outdir + '/' + simstr + '_pc_' + str(pcid) + '.npz'
  print 'host ' , pcid, ' saving to ' , fn
  ne,ni,szfast = 0,0,0
  lE,lI=[],[]
  for c in ce:
    if ice(c.ty):
      lI.append(c.ID)
      ni += 1
    else:
      lE.append(c.ID)
      ne += 1
    szfast = int(c.soma_volt.size())
  lE=array(lE) # lE is list of E cell IDs from this host
  lI=array(lI) # Li is list of I cell IDs from this host
  soma_volt = zeros((ne,szfast)); Adend3_volt = zeros((ne,szfast)); Bdend_volt=zeros((ne,szfast));
  soma_voltI = zeros((ni,szfast));
  cdx = 0; idx = 0;
  for c in ce:
    if ice(c.ty):
      soma_voltI[idx,:] = c.soma_volt.to_python()
      idx += 1
      continue
    soma_volt[cdx,:] = c.soma_volt.to_python()
    Adend3_volt[cdx,:] = c.Adend3_volt.to_python()
    Bdend_volt[cdx,:] = c.Bdend_volt.to_python()
    cdx += 1
  numpy.savez(fn,lctyID=array(lctyID),lX=array(lX),lY=array(lY),lZ=array(lZ),vt=vt.as_numpy(),lspks=lspks,lids=lids,lE=lE,lI=lI,Adend3_volt=Adend3_volt,Bdend_volt=Bdend_volt)

pc.barrier()

####################################################################################

if saveout: # save the sim data
  if pcID == 0: print 'saving data'
  savedata(simstr,pcID)

if saveconns and pcID == 0: # save connectivity data
  print 'saving connections'
  h.batch_flag = 1 # in case file already exists
  connsnq.tog("DB")
  connsnq.sv(outdir + '/' + simstr + 'connsnq.nqs')

pc.runworker()
pc.done()

if nhosts > 1: h.quit() # this means was likely running in batch mode
