from pfac.fac import *
import sys, os
import numpy as np

z = int(sys.argv[1])
k = int(sys.argv[2])

if not os.path.exists('dav'):
    os.system('mkdir dav')
a = ATOMICSYMBOL[z]
p = 'dav/%s%02d'%(a,k)

SetUTA(1)

SetAtom(a)

SetOption('structure:full_name', 1)
Closed('1s1 2*1 3*1 4[s,p,d]1 5[s,p]1')
v = k-54

gs = []
for n5d in range(3):
    for n5f in range(2):
        for n6s in range(3):
            for n6p in range(2):
                for n6d in range(2):
                    n4f = v - (n5d+n5f+n6s+n6p+n6d)
                    if n4f < 0:
                        continue
                    nn = (n4f,n5d,n5f,n6s,n6p,n6d)
                    gs.append('4f%d.5d%d.5f%d.6s%d.6p%d.6d%d'%nn)
                    Config(gs[-1],'4f%d 5d%d 5f%d 6s%d 6p%d 6d%d'%nn)

ConfigEnergy(0)
OptimizeRadial(gs)
ConfigEnergy(1)
GetPotential(p+'a.pot')

Config(1, 'ge', gs, '4f1 5d1 5f1 6s1 6p1 6d1', 4, 6, 0, 5)
ListConfig(p+'a.cfg')

gs.append('ge')
Structure(p+'b.en', gs)

MemENTable(p+'b.en')
PrintTable(p+'b.en', p+'a.en')
