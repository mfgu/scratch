from pfac.fac import *
import sys

md = int(sys.argv[1])
ConvertToSFAC('m.sf')
InitializeMPI(2)
a = 'Fe'
p = '%s10'%a
nmin=3
nmax=10
smax=4
SetAtom(a)
gv=['g2']
Config('g2', '1s2 2*8')
for n in range(3,smax+1):
    gn = 'g%d'%n
    gv.append(gn)
    Config(gn, '1s2 2*7 %d*1'%n)

WallTime('opt')
ConfigEnergy(0)
OptimizeRadial('g2')
ConfigEnergy(1)

iv=[]
for n in range(3, nmax+1):    
    gn = 'i%d'%n
    iv.append(gn)
    Config(gn, '1s2 2*8 %d*1'%n)
dv=[]
nv=[]
mv=[]
for n in [3,4,5]:
    for m in range(max(nmin,n),nmax+1):
        gn = 'd%d.%d'%(n,m)
        dv.append(gn)
        nv.append(n)
        mv.append(m)
        if m == n:
            Config(gn, '1s2 2*7 %d*2'%m)
        else:
            Config(gn, '1s2 2*7 %d*1 %d*1'%(n,m))

WallTime('en')
for g in gv+iv+dv:
    Structure(p+'b.en', [g])

MemENTable(p+'b.en')
PrintTable(p+'b.en', p+'a.en')

WallTime('tr')
for m in [-1, 1, -2, 2]:
    TRTable(p+'b.tr', ['g2'], ['g3'], m)
    TRTable(p+'b.tr', ['g3'], ['g3'], m)

m = -1
for i0 in range(4,smax+1):
    TRTable(p+'b.tr', ['g2'], [gv[i0-2]], m)
    TRTable(p+'b.tr', ['g3'], [gv[i0-2]], m)
    for i1 in range(4,i0+1):
        TRTable(p+'b.tr', [gv[i1-2]], [gv[i0-2]], m)
        
for i0 in range(len(iv)):
    for i1 in range(i0, len(iv)):
        TRTable(p+'b.tr', [iv[i0]], [iv[i1]], m)
for i in range(len(dv)):
    TRTable(p+'b.tr', [iv[nv[i]-3]], [dv[i]], m)
    if mv[i] != nv[i]:
        TRTable(p+'b.tr', [iv[mv[i]-3]], [dv[i]], m)
    for j in range(len(dv)):
        if mv[j] == mv[i] and nv[j] < nv[i]:
            TRTable(p+'b.tr', [dv[j]], [dv[i]], m)
PrintTable(p+'b.tr', p+'a.tr')

WallTime('ce')
for i0 in range(3,smax+1):
    if md > 0:
        CETableMSub(p+'b.ce', ['g2'], [gv[i0-2]])
    else:
        CETable(p+'b.ce', ['g2'], [gv[i0-2]])
PrintTable(p+'b.ce', p+'a.ce')

WallTime('ai')
for i in range(len(dv)):
    if md > 0:
        AITableMSub(p+'b.ai', [dv[i]], ['g2'])
        AITableMSub(p+'b.ai', [dv[i]], ['g3'])
        if nv[i] > 3:
            AITableMSub(p+'b.ai', [dv[i]], ['g4'])
    else:
        AITable(p+'b.ai', [dv[i]], ['g2'])
        AITable(p+'b.ai', [dv[i]], ['g3'])
        if nv[i] > 3:
            AITable(p+'b.ai', [dv[i]], ['g4'])
PrintTable(p+'b.ai', p+'a.ai')

WallTime('done')
FinalizeMPI()
CloseSFAC()
