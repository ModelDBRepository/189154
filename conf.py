import ConfigParser
import io

# default config as string
def_config = """
[seed]
iseed = 1234
wseed = 4321
pseed = 4321
[netsyn]
NMAMREE = 0.1
NMAMREI = 0.1
mGLURR = 7.5
GB2R = 7.5
rdmsec = 1
nmfracca = 0.13
[chan]
ihginc = 2.0
iark2fctr = 1.0
iark4 = 0.008
erevh = -30.0
h_lambda = 325.0
h_gbar = 0.0025
fs_h_gbar = 0.00002
lts_h_gbar = 0.15
cagk_gbar = 0.0001
ikc_gkbar = 0.003
nax_gbar = 0.081
kdr_gbar = 0.021
kap_gbar = 0.3
kdmc_gbar = 0.00085
km_gmax = 0.1
cabar = 0.005
lts_cabar = 1.0
[cada]
taur = 5
[run]
indir = data
outdir = data
tstop = 2000.0
dt = 0.1
saveout = 1
simstr = 15dec29_B
statestr = 15apr20_net_S3
dorun = 1
doquit = 0
dodraw = 0
verbose = 0
recdt = 10.0
recvdt = 1.0
binsz = 5
saveconns = 0
[rxd]
CB_frate=5.5
CB_brate=0.0026
CB_init=0.2
gip3 = 120400.0
gserca = 4.0
gleak = 3.0
cacytinit = 100e-6
caerinit = 1.25
caexinit = 0.0
spaceum = 0.0
nsubseg = 0
subsegum = 0.0
v1ryr = 100.0
[net]
scale=1.0
IIGain = 0.1
IEGain = 0.15
EIGainFS = 0.15
EIGainLTS = 0.15
EEGain = 0.25
[stim]
EXGain = 1.0
noise = 1
ip3_stim = 0.0
ip3_stimT = 10000.0
sgrhzNMI = 600.0
sgrhzNME = 300.0
sgrhzEE = 800.0
sgrhzEI = 1600.0
sgrhzIE = 150.0
sgrhzII = 150.0
sgrhzMGLURE = 0.0
sgrhzGB2 = 0.0
"""

# write config file starting with defaults and new entries
# specified in section (sec) , option (opt), and value (val)
# saves to output filepath fn
def writeconf (fn,sec,opt,val):
  conf = ConfigParser.ConfigParser()
  conf.optionxform = str
  conf.readfp(io.BytesIO(def_config)) # start with defaults
  # then change entries by user-specs
  for i in xrange(len(sec)): conf.set(sec[i],opt[i],val[i])
  # write config file
  with open(fn, 'wb') as cfile: conf.write(cfile)

# read config file
def readconf (fn="physiol.cfg"):

  config = ConfigParser.ConfigParser()
  config.optionxform = str
  config.read(fn)

  def conffloat (base,var,defa): # defa is default value
    val = defa
    try: val=config.getfloat(base,var)
    except: pass
    return val

  def confint (base,var,defa):
    val = defa
    try: val=config.getint(base,var)
    except: pass
    return val

  def confstr (base,var,defa):
    val = defa
    try: val = config.get(base,var)
    except: pass
    return val

  d = {}

  d['iseed'] = confint("seed","iseed",1234)
  d['wseed'] = confint("seed","wseed",4321)
  d['pseed'] = confint("seed","pseed",4321)
  d['NMAMREE'] = conffloat("netsyn","NMAMREE",0.1)
  d['NMAMREI'] = conffloat("netsyn","NMAMREI",0.1)
  d['mGLURR'] = conffloat("netsyn","mGLURR",7.5)
  d['GB2R'] = conffloat("netsyn","GB2R",7.5)
  d['nmfracca'] = conffloat("netsyn","nmfracca", 0.13)
  d['rdmsec'] = confint("netsyn","rdmsec", 1)
  d['erevh'] = conffloat("chan","erevh",-30.0)
  d['h_lambda'] = conffloat("chan","h_lambda",325.0)
  d['h_gbar'] = conffloat("chan","h_gbar",0.0025)
  d['fs_h_gbar'] = conffloat("chan","fs_h_gbar",0.00002)
  d['lts_h_gbar'] = conffloat("chan","lts_h_gbar",0.15)
  d['cagk_gbar'] = conffloat("chan","cagk_gbar",0.0001)
  d['ikc_gkbar'] = conffloat("chan","ikc_gkbar",0.003)
  d['nax_gbar'] = conffloat("chan","nax_gbar",0.081)
  d['kdr_gbar'] = conffloat("chan","kdr_gbar",0.021)
  d['kap_gbar'] = conffloat("chan","kap_gbar",0.3)
  d['kdmc_gbar'] = conffloat("chan","kdmc_gbar",0.00085)
  d['km_gmax'] = conffloat("chan","km_gmax",0.1)
  d['ihginc'] = conffloat("chan","ihginc", 2.0)
  d['iark2fctr'] = conffloat("chan", "iark2fctr",1.0)
  d['iark4'] = conffloat("chan", "iark4",0.008)
  d['cabar'] = conffloat("chan","cabar",0.005)
  d['lts_cabar'] = conffloat("chan","lts_cabar",1.0)
  d['taurcada'] = conffloat("cada", "taur", 5.0)
  d['outdir'] = confstr("run","outdir", "data")
  d['indir'] = confstr("run","indir", "data")
  d['tstop'] = conffloat("run","tstop", 2000.0)
  d['dt'] = conffloat("run","dt",0.1)
  d['saveout'] = conffloat("run","saveout",1)
  d['simstr'] = confstr("run","simstr","15dec29_B")
  d['statestr'] = confstr("run","statestr","15apr20_net_S3")
  d['dorun'] = confint("run","dorun",1)
  d['recdt'] = conffloat("run","recdt",10.0)
  d['recvdt'] = conffloat("run","recvdt",1.0)
  d['binsz'] = conffloat("run","binsz",5)

  for k in ['saveconns','doquit','verbose','dodraw']: d[k] = confint("run",k,0)

  d['CB_frate'] = conffloat("rxd","CB_frate", 5.5)
  d['CB_brate'] = conffloat("rxd","CB_brate", 0.0026)
  d['CB_init'] = conffloat("rxd","CB_init", 0.2)
  d['gip3'] = conffloat("rxd","gip3",120400.0)
  d['gserca'] = conffloat("rxd","gserca",4.0)
  d['gleak'] = conffloat("rxd","gleak",3.0)
  d['caerinit'] = conffloat("rxd","caerinit",1.25)
  d['cacytinit'] = conffloat("rxd","cacytinit",100e-6)
  d['caexinit'] = conffloat("rxd","caexinit",0.0)
  d['spaceum'] = conffloat("rxd","spaceum",0.0)
  d['nsubseg'] = confint("rxd","nsubseg",0)
  d['subsegum'] = conffloat("rxd","subsegum",0.0)
  d['v1ryr'] = conffloat("rxd","v1ryr",100.0)
  d['scale'] = conffloat("net","scale",1.0)
  d['IIGain'] = conffloat("net","IIGain",0.1)
  d['IEGain'] = conffloat("net","IEGain",0.15)
  d['EIGainFS'] = conffloat("net","EIGainFS",0.15)
  d['EIGainLTS'] = conffloat("net","EIGainLTS",0.15)
  d['EEGain'] = conffloat("net","EEGain",0.25)
  d['EXGain'] = conffloat("stim","EXGain",1.0)
  d['noise'] = conffloat("stim","noise",1.0)
  d['ip3_stim'] = conffloat("stim","ip3_stim",0.0)
  d['ip3_stimT'] = conffloat("stim","ip3_stimT",10000.0)
  d['sgrhzNME'] = conffloat("stim","sgrhzNME",300.0)
  d['sgrhzNMI'] = conffloat("stim","sgrhzNMI",600.0)
  d['sgrhzEE'] = conffloat("stim","sgrhzEE",800.0)
  d['sgrhzIE'] = conffloat("stim","sgrhzIE",150.0)
  d['sgrhzEI'] = conffloat("stim","sgrhzEI",1600.0)
  d['sgrhzII'] = conffloat("stim","sgrhzII",150.0)
  d['sgrhzMGLURE'] = conffloat("stim","sgrhzMGLURE",0.0)
  d['sgrhzGB2'] = conffloat("stim","sgrhzGB2",0.0)

  return d

