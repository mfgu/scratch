from pfac.pol import *
import sys
from multiprocessing import Pool

#ConvertToSPOL('p.sf')

m = int(sys.argv[1])
np = int(sys.argv[2])
a = 'Fe'
p = '%s10'%a
e0 = 300.0
nn = 1000
sig = 1.7
de = sig/2.0
etr = 100.0
e = [e0 + i*de for i in range(nn)]

def ploop(i0):
    SetMaxLevels(89)
    SetMLevels(p+'b.en', p+'b.tr')
    for i in range(i0,nn,np):
        Print('%d %10.3E'%(i, e[i]))
        SetMLevels()
        SetEnergy(e[i], sig)
        if m != 2:
            SetMCERates(p+'b.ce')
        if m != 1:
            SetMAIRates(p+'b.ai')
        SetDensity(1.0)
        PopulationTable('pop%de%03d.txt'%(m,i))
        Orientation(etr)
        PolarizationTable('pol%de%03d.txt'%(m,i), '', 700, 850, 1e-8)

p = Pool(processes=np)
p.map(ploop, range(np))

#CloseSPOL()
