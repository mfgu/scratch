from pfac.fac import *
import sys, os
import numpy as np
from pfac import rfac

z = int(sys.argv[1])
k = int(sys.argv[2])
e0 = float(sys.argv[3])
e1 = float(sys.argv[4])
e2 = float(sys.argv[5])
if (len(sys.argv)>6):
    csf = sys.argv[6]

if not csf is None:
    ConvertToSFAC(csf)
    
InitializeMPI(8)

if not os.path.exists('dci'):
    os.system('mkdir dci')
a = ATOMICSYMBOL[z]
p = 'dci/%s%02d'%(a,k)

SetAtom(a)

SetOption('structure:full_name', 1)
Closed('1s1 2*1 3*1 4[s,p,d]1 5[s,p]1')

r = rfac.FLEV('dav/%s%02da.en'%(a,k))
w0 = np.where((r.e-r.e0)<e0)[0]
w1 = np.where((r.e-r.e0)<e1)[0]
cs0 = np.unique(r.s[w0])
cs = np.unique(r.s[w1])

gs = []
gs0 = []
gs1 = []
cs1 = []
for i in range(len(cs)):
    gs.append('g%02d'%i)
    if (cs[i] in cs0):
        gs0.append(gs[-1])
    else:
        gs1.append(gs[-1])
        cs1.append(cs[i])
    Config(gs[-1], cs[i].replace('.', ' '))

for i in range(len(gs0)):
    print('gs0: %s %s'%(gs0[i],cs0[i]))
for i in range(len(gs1)):
    print('gs1: %s %s'%(gs1[i],cs1[i]))
    
ListConfig(p+'a.cfg')
WallTime('opt')
OptimizeRadial(gs0)
GetPotential(p+'a.pot')

WallTime('en cs0')
Structure(p+'b.en', gs0)
for i in range(len(gs1)):
    WallTime('en '+cs1[i])
    Structure(p+'b.en', [gs1[i]])

MemENTable(p+'b.en')
PrintTable(p+'b.en', p+'a.en')

SetOption('transition:lower_emax', e2)
SetOption('transition:tr_emin', 0.01)
SetOption('transition:tr_emax', 25.0)
WallTime('tr')
TRTable(p+'b.tr', gs, gs, -1)
PrintTable(p+'b.tr', p+'a.tr')

FinalizeMPI()

if not csf is None:
    CloseSFAC()
