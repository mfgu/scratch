from pfac.pol import *
import sys

#ConvertToSPOL('p.sf')

m = int(sys.argv[1])

a = 'Fe'
p = '%s10'%a
e0 = 300.0
nn = 1000
sig = 1.6
de = sig/2.0
e = [e0 + i*de for i in range(nn)]
SetMaxLevels(89)
SetMLevels(p+'b.en', p+'b.tr')
for i in range(0,nn):
    Print('%d %10.3E'%(i, e[i]))
    SetMLevels()
    SetEnergy(e[i], sig)
    if m != 2:
        SetMCERates(p+'b.ce')
    if m != 1:
        SetMAIRates(p+'b.ai')
    SetDensity(1.0)
    PopulationTable('pop%de%03d.txt'%(m,i))
    Orientation()
    PolarizationTable('pol%de%03d.txt'%(m,i))

#CloseSPOL()
