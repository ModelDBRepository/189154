import sys
print sys.path
import os
import string
from neuron import *
from datetime import datetime
h("strdef simname, allfiles, simfiles, output_file, datestr, uname, osname, comment")
h.simname=simname = "dystdemo"
h.allfiles=allfiles = "pyinit.py geom.py mpisim.py"
h.simfiles=simfiles = "pyinit.py geom.py mpisim.py"
h("runnum=1")
runnum = 1.0
h.datestr=datestr = "16jun21"
h.output_file=output_file = "data/16jun21.0"
h.uname=uname = "x86_64"
h.osname=osname="linux"
h("templates_loaded=0")
templates_loaded=0
h("xwindows=1.0")
xwindows = 1.0
h.xopen("nrnoc.hoc")
h.xopen("init.hoc")
CTYPi = 60.0; STYPi = 18.0
from pyinit import *
from neuron import h, gui
from vector import *
from nqs import *
from labels import *
import random
from pylab import *
import ConfigParser

tl = tight_layout
ion()

rcParams['lines.markersize'] = 15
rcParams['lines.linewidth'] = 4
rcParams['font.size'] = 25

defbinsz=10
stimdel = 500 # delay from stim to measured firing rate peak response

defnCPU=8
print 'Using ' , defnCPU, ' cpus by default.'

config = ConfigParser.ConfigParser()

# determine config file name
def setfcfg ():
  fcfg = "physiol.cfg" # default config file name
  for i in xrange(len(sys.argv)):
    if sys.argv[i].endswith(".cfg") and os.path.exists(sys.argv[i]):
      fcfg = sys.argv[i]
  print "config file is " , fcfg
  return fcfg

fcfg=setfcfg() # config file name

config.read(fcfg)

def conffloat (base,var): return float(config.get(base,var))
def confint (base,var): return int(config.get(base,var))
def confstr (base,var): return config.get(base,var)

tstart = 0
tstop = conffloat('run','tstop') + tstart
binsz = conffloat("run","binsz") # bin size (in milliseconds) for MUA
recdt = conffloat('run','recdt')
recvdt = conffloat('run','recvdt')
vtslow=Vector(); vtslow.indgen(tstart,tstop-recdt,recdt); vtslow=numpy.array(vtslow.to_python())
vtfast=Vector(); vtfast.indgen(tstart,tstop-recvdt,recvdt); vtfast=numpy.array(vtfast.to_python())
simstr = config.get('run','simstr')

def makefname (simstr,pcid): return 'data/' + simstr + '_pc_' + str(pcid) + '.npz'

ix,ixe=[1e9 for i in xrange(CTYPi)],[-1e9 for i in xrange(CTYPi)]
numc=[0 for i in xrange(CTYPi)]
allcells,icells,ecells=0,0,0

# reload the variables from config file
def reloadvars ():
  global tstart,loadstate,tstop,binsz,recdt,recvdt,vtslow,vtfast,simstr,numc,allcells,ecells,icells,startt
  tstart = conffloat('run','loadtstop')
  loadstate = confint('run','loadstate')
  if loadstate == 0: tstart = 0 # only use previous end time if loading state
  tstop = conffloat('run','tstop') + tstart
  binsz = conffloat("run","binsz") # bin size (in milliseconds) for MUA
  recdt = conffloat('run','recdt')
  recvdt = conffloat('run','recvdt')
  vtslow=Vector(); vtslow.indgen(tstart,tstop-recdt,recdt); vtslow=numpy.array(vtslow.to_python())
  vtfast=Vector(); vtfast.indgen(tstart,tstop-recvdt,recvdt); vtfast=numpy.array(vtfast.to_python())
  simstr = config.get('run','simstr')
  numc=[0 for i in xrange(CTYPi)]
  allcells,ecells,icells=0,0,0

# setup start/end indices for cell types
def makeix (lctyID):
  global ix,ixe,allcells,ecells,icells
  allcells,ecells,icells=0,0,0
  for i in xrange(CTYPi):
    ix[i]=1e9
    ixe[i]=-1e9
    numc[i]=0
  for i in xrange(len(lctyID)):
    ty = lctyID[i]
    numc[ty]+=1
    allcells+=1
    ix[ty] = min(ix[ty],i)
    ixe[ty] = max(ixe[ty],i)
    if h.ice(ty): icells+=1
    else: ecells+=1

lfastvar = ['volt', 'iAM', 'iNM', 'iGB', 'iGA','ina','ik','ica','ih']

# reads data from files saved by mpisim.py - must get nhost correct
def readdat (simstr,nhost):
  ld = {} # dict for data from a single host
  for pcid in xrange(nhost):
    fn = makefname(simstr,pcid)
    ld[pcid]=numpy.load(fn)
  return ld
  
# load/concat data from mpisim.py (which saves data from each host separately) - must get nhost correct
def loaddat (simstr,nhost,quiet=False):
  ld = readdat(simstr,nhost) # dict for data from a single host
  lids,lspks=array([]),array([]) # stitch together spike times
  for pcid in xrange(nhost):
    lids=concatenate((lids,ld[pcid]['lids']))
    lspks=concatenate((lspks,ld[pcid]['lspks']))
  ldout={} # concatenated output from diff hosts
  ldout['lids']=lids; ldout['lspks']=lspks;
  for k in ['lctyID','vt','lX', 'lY', 'lZ']: ldout[k]=ld[0][k][:]
  makeix(ldout['lctyID'])
  ncells = len(ld[0]['lctyID'])
  vt = ld[0]['vt'][:]
  locvarlist = []
  for loc in ['soma','Adend3']:
    for var in ['volt']:
      locvarlist.append(loc + '_' + var)
  locvarlist.append('Bdend_volt')
  for k in locvarlist:
    if not k in ld[0]: continue # skip?
    if not quiet: print k
    if k.endswith('volt') or lfastvar.count(k.split('_')[1])>0: sz = len(vtfast)
    else: sz = len(vtslow)
    ldout[k] = zeros((ncells,sz))
    ldk = ldout[k]
    for pcid in xrange(nhost):
      d = ld[pcid][k]
      lE = ld[pcid]['lE']; 
      for row,cdx in enumerate(lE): ldk[cdx,:] = d[row,0:sz]
      if k == 'soma_volt': # special case for soma voltage of I cells
        lI=ld[pcid]['lI']; dI=ld[pcid]['soma_voltI']; row=0
        for row,cdx in enumerate(lI): ldk[cdx,:] = dI[row,0:sz]
  for pcid in xrange(nhost):
    del ld[pcid].f
    ld[pcid].close()
    del ld[pcid]
  return ldout

# colors for spikes for raster -- based on cell type
def getcolors (lspks,lids,lctyID):
  lclrs = []
  for i in xrange(len(lids)):
    ty = lctyID[int(lids[i])]
    if IsLTS(ty):
      lclrs.append('b')
    elif ice(ty):
      lclrs.append('g')
    else:
      lclrs.append('r')
  return lclrs
			
#
def drawraster (ld,sz=2):
  lspks,lids,lctyID=ld['lspks']/1e3,ld['lids'],ld['lctyID']
  try:
    lclrs = ld['lclrs']
  except:    
    lclrs = getcolors(lspks,lids,lctyID)
    ld['lclrs'] = lclrs
  scatter(lspks,lids,s=sz**2,c=lclrs,marker='+')
  allcells = len(lctyID)
  xlim((tstart/1e3,tstop/1e3)); ylim((0,allcells)); tight_layout(); xlabel('Time (s)',fontsize=30); 
  ax = gca()
  ax.set_yticks([])

# get cells in vid if they're of specified type (in lty)
def getall (ix,ixe,lty=[E2,E4,E5R,E5B,E5P,E6,E6C]):
  vact = Vector()
  for ct in lty:
    vtmp = Vector(); vtmp.indgen(ix[ct],ixe[ct],1); vact.append(vtmp)
  return vact

# get an array of times of interest
def gettimes (tstart,tstop,binsz):
  vt = h.Vector()
  vt.indgen(tstart,tstop,binsz)
  return vt

# get an NQS with spikes (snq)
def getsnq (ld):
  lspks,lids,lctyID=ld['lspks'],ld['lids'],ld['lctyID']
  snq = NQS('id','t','ty','ice')
  snq.v[0].from_python(lids)
  snq.v[1].from_python(lspks)
  for i in xrange(len(lids)):
    snq.v[2].append(lctyID[int(lids[i])])
    snq.v[3].append(h.ice(lctyID[int(lids[i])]))
  ld['snq']=snq
  return snq

# print spike rates in each population in the periods of interest (baseline, signal, recall)
def prspkcount (snq,tstart,tstop,binsz,times=None,quiet=False):
  snq.verbose=0
  if times is None:
    times = gettimes(tstart,tstop,binsz) # periods of interest
  lfa=[]
  vact = getall(ix,ixe,lty=[E2,E4,E5R,E5B,E5P,E6,E6C])
  for startt in times:
    endt = startt + binsz
    if endt > tstop: break
    fa = round(1e3*snq.select("id","EQW",vact,"t","[]",startt,endt) / ((endt-startt)*vact.size()),2)
    lfa.append(fa)
    if not quiet:
      stro = "t=" + str(startt) + "-" + str(endt) + ". E:" + str(fa) + " Hz."
      for ct in [I2, I2L, I5, I5L, I6, I6L]:
        fb = round(1e3*snq.select("id","[]",ix[ct],ixe[ct],"t","[]",startt,endt) / ((endt-startt)*numc[ct]),2)
        stro += ' ' + CTYP[ct] + ':' + str(fb)
      fi = round(1e3*snq.select("ice",1,"t","[]",startt,endt) / ((endt-startt)*icells),2)
      stro += ' I:' + str(fi)
      for ct in [E2, E4, E5R, E5B, E5P, E6, E6C]:
        fb = round(1e3*snq.select("id","[]",ix[ct],ixe[ct],"t","[]",startt,endt) / ((endt-startt)*numc[ct]),2)
        stro += ' ' + CTYP[ct] + ':' + str(fb)
      fe = round(1e3*snq.select("ice",0,"t","[]",startt,endt) / ((endt-startt)*ecells),2)
      stro += ' E:' + str(fe)
      print stro
  snq.verbose=1
  return lfa

# getfnq - make an NQS with ids, firing rates, types
def getfnq(snq,lctyID,skipms=200):
  snq.verbose=0; snq.tog("DB"); 
  fnq = h.NQS("id","freq","ty")
  tf = tstop - skipms # duration we're considering for frequency calc
  for i in xrange(allcells):
    n = float( snq.select("t",">",skipms,"id",i) ) # number of spikes
    fnq.append(i, n*1e3/tf, lctyID[i])
  snq.verbose=1
  return fnq

# pravgrates - print average firing rates using fnq
def pravgrates(snq,skipms=500):
  snq.verbose=0
  for ty in xrange(CTYPi):
    if numc[ty] < 1: continue
    print CTYP[ty], ' avg rate = ', getrate(snq,ty,500.0)[2]
  snq.verbose=1

###########################################################3

# reload variables, reset config file name
def myreload (Simstr,ncpu=defnCPU, quiet=False):
  global fcfg,simstr
  if not quiet: print 'simstr was ' , simstr
  simstr = Simstr
  print 'simstr: ' , simstr
  fcfg = 'backupcfg/' + simstr + '.cfg'
  if len(config.read(fcfg))==0:
    print 'myreload ERRA: could not read ', fcfg, '!!!'
    return None
  reloadvars()
  if not quiet: print fcfg
  try:
    ld=loaddat(simstr,ncpu,quiet);
    return ld
  except:
    print 'could not reload data from' , Simstr
    return None

ld=None
  
#
def mydrawrast (lld):
  for i,ld in enumerate(lld):
    figure(); 
    drawraster(ld,sz=4); title(lsimstr[i])
    xlabel('')
    xlim((0.5,2))  
    fsz= 20
    ax = gca(); 
    ax.set_yticks([]); ax.set_xticks([]); ax.set_xlabel('');
    tx=-0.03
    lty = [(ix[ct]+(ixe[ct]-ix[ct])*0.3)/allcells for ct in [E2,E4,E5R,E5B,E5P,E6C,E6]]
    for ty,lbl in zip(lty,['E2','E4','E5a','E5b','E5P','E6C','E6']): addtext(1,1,[1],[lbl],tx=tx,ty=ty,c='r',fsz=fsz+5)
    lty = [(ix[ct]+(ixe[ct]-ix[ct])*0.075)/allcells for ct in [I2,I5,I6]]
    for ty,lbl in zip(lty,['I2','I5','I6']): addtext(1,1,[1],[lbl],tx=tx,ty=ty,c='g',fsz=fsz)
    lty = [(ix[ct]+(ixe[ct]-ix[ct])*0.25)/allcells for ct in [I2L,I5L,I6L]]
    for ty,lbl in zip(lty,['I2L','I5L','I6L']): addtext(1,1,[1],[lbl],tx=tx,ty=ty,c='b',fsz=fsz)

#
def addtext (row,col,lgn,ltxt,tx=-0.025,ty=1.03,c='k',fsz=30):
  for gn,txt in zip(lgn,ltxt):
    if row == col and col == 1:
      ax = gca()
    else:
      ax = subplot(row,col,gn)
    text(tx,ty,txt,fontweight='bold',transform=ax.transAxes,fontsize=fsz,color=c);

#
def Gaddtext (G,row,col,txt,tx=-0.025,ty=1.03,c='k',fsz=30):
  ax = subplot(G[row,col])
  text(tx,ty,txt,fontweight='bold',transform=ax.transAxes,fontsize=fsz,color=c);

def naxbin (ax,nb): ax.locator_params(nbins=nb);

lsimstr = ['d2', 'dystonia', 'latchup']
lld,lsnq = [],[]
for simstr in lsimstr:
  ld = loaddat(simstr,defnCPU)
  lld.append(ld)
  snq = getsnq(ld)
  lsnq.append(snq)

mydrawrast(lld)

