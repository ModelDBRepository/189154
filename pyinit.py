# $Id: pyinit.py,v 1.2 2012/02/15 16:22:52 samn Exp $ 

from neuron import h
import os
import sys
import datetime
import shutil
import pickle
from math import sqrt, pi
import numpy
import types

h("objref p")
h("p = new PythonObject()")

try:
    import pylab
    from pylab import plot, arange, figure
    my_pylab_loaded = True
except ImportError:
    print "Pylab not imported"
    my_pylab_loaded = False

def htype (obj): st=obj.hname(); sv=st.split('['); return sv[0]
def secname (obj): obj.push(); print h.secname() ; h.pop_section()
def psection (obj): obj.push(); print h.psection() ; h.pop_section()

allsecs=None #global list containing all NEURON sections, initialized via mkallsecs

# still need to generate a full allsecs
def mkallsecs ():
  """ mkallsecs - make the global allsecs variable, containing
      all the NEURON sections.
  """
  global allsecs
  allsecs=h.SectionList() # no .clear() command
  roots=h.SectionList()
  roots.allroots()
  for s in roots:
    s.push()
    allsecs.wholetree()
  return allsecs

#forall syntax - c gets executed, allsecs has Sections
def forall (c):
    """ NEURON forall syntax - iterates through all the sections available
        note that there's a dummy loop variable called s used in this function,
        so any command that needs to access a section should be via s.
        example: forall('print s.name()') , will print all the section names.
        Also note that this function uses a global list, 'allsecs', which may
        need to get re-initialized when new sections are created, via the mkallsecs
        function above.
    """
    global allsecs
    if (type(allsecs)==types.NoneType): mkallsecs()
    for s in allsecs: exec(c)

#forsec syntax - executes command for each section who's name
# contains secname as a substring
def forsec (secref="soma",command=""): 
    """ NEURON forsec syntax - iterates over all sections which have a substring
        in their names matching secref argument. command is executed if match found.
        this function also utilizes the allsecs global variable.
    """
    global allsecs
    if (type(allsecs)==types.NoneType): mkallsecs()
    if (type(secref)==types.StringTypes[0]):
        for s in allsecs:
            if s.name().count(secref) > 0:
                exec(command)
    else:
        for s in secref: exec(command)
