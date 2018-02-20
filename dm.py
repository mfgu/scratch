from optparse import OptionParser
from pfac.fac import *
from time import time

t0 = time()

p = OptionParser()
p.add_option('-z', dest='z', type='int',
             default=26, help='atomic number')
p.add_option('-n', dest='n', type='int',
             default=1, help='number of M-shell electrons')
p.add_option('-p', dest='np', type='int',
             default=1, help='number of processors')
p.add_option('-m', dest='nmax', type='int',
             default=5, help='max principle quantum number')
p.add_option('-i', dest='ni', type='int',
             default=1, help='inner shell excitation')
p.add_option('-r', dest='nr', type='int',
             default=0, help='job id for the n1 split run')
p.add_option('-s', dest='nsp', type='int',
             default=4, help='number of n1s in split jobs')
p.add_option('-t', dest='ntr', type='int',
             default=3, help='mbpt transition rate option')
p.add_option('-k', dest='nk', type='int',
             default=0, help='k shell max excitation')
p.add_option('-c', dest='csf', type='int',
             default=0, help='convert to sfac input file')
p.add_option('-a', dest='mm', type='float',
             default=0, help='maximum memory usage for radial integra caching')
p.add_option('-d', dest='dry', type='int',
             default=0, help='dry run, stop after configuration setup')

(opts, args) = p.parse_args()

print opts

i = opts.nr
asym = ATOMICSYMBOL[opts.z]

if opts.csf > 0:
    if i >= 0:
        ConvertToSFAC('dm%d.sf'%i)
    else:
        ConvertToSFAC('dm.sf')
    
InitializeMPI(opts.np)
SetAtom(asym)

n = opts.n
nmax = opts.nmax
pref='%s%02d'%(asym, n+10)
if i >= 0:
    p0 = '%si%d'%(pref,i)
else:
    p0 = pref
    
gc=[]
gc.append('g3')
Config('g3', '1s2 2*8 3*%d'%n)
for m in range(4,nmax+1):
    gc.append('g%d'%m)
    Config('g%d'%m, '1s2 2*8 3*%d %d*1'%(n-1, m))
if opts.ni > 0:
    gc.append('i3')
    Config('i3', '1s2 2*7 3*%d'%(n+1))
    for m in range(4, nmax+1):
        gc.append('i%d'%m)
        Config('i%d'%m, '1s2 2*7 3*%d %d*1'%(n, m))

ListConfig(p0+'a.cfg')

n1 = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20, 22, 24, 28, 32, 38, 43, 50, 65, 80, 100, 125, 150, 175, 200, 225, 250]
nn = len(n1)
if (opts.nsp <= 0):
    opts.nsp = nn
ns = range(0, nn, opts.nsp)
ns.append(nn)

n2 = range(8)+[8, 13, 20, 28, 38, 50, 65, 80, 100]

eps = 1e-4

SetVP(103)
SetMS(3, 3)
SetSE(-1, 61)
SetBreit(-1, 1, 5)

PrintNucleus()
PrintQED()
Print('ns=%d'%(len(ns)-1))
Print(ns)
if i >= 0:
    print n1[ns[i]:ns[i+1]]
    print n2
    
if opts.dry > 0:
    exit(0)
    
if (opts.mm > 0):
    LimitArray(-1, opts.mm)
    
iop=0
if i >= 0:
    OptimizeRadial(gc[iop])
    SetBoundary(nmax, eps, 1e30)
    ReinitRadial(0)
    SetRadialGrid(3000, 1.1, -1e30, 0.0)
    
SetPotentialMode(10)
OptimizeRadial(gc[iop])
SetBoundary(nmax, eps, 1e30)
GetPotential(p0+'a.pot')

TransitionMBPT(opts.ntr, opts.ntr)
if i >= 0:
    StructureMBPT(0, 0, eps)
    if opts.nsp >= nn:
        StructureMBPT(1)
    if opts.nk >= 0:
        StructureMBPT('1s', opts.nk)
    StructureMBPT(p0+'b.en', p0+'b.ham', gc,
                  [-1,-1,-1,-1,-1,-1,-1,-1,15,15,15,15],
                  n1[ns[i]:ns[i+1]], n2, len(gc))
    MemENTable(p0+'b.en')
    PrintTable(p0+'b.en', p0+'a.en')
else:
    h = [pref+'i%02db.ham'%x for x in range(len(ns)-1)]
    StructureMBPT(1)
    if (opts.ntr > 0):
        TransitionMBPT(p0+'b.tr', gc, gc)
    StructureMBPT(p0+'b.en', pref+'a.ham', h, gc, len(gc))
    MemENTable(p0+'b.en')
    PrintTable(p0+'b.en', p0+'a.en')
    if (opts.ntr > 0):
        PrintTable(p0+'b.tr', p0+'a.tr')

t1 = time()
Print('all done %d %10.3E'%(i,t1-t0))
FinalizeMPI()
if opts.csf > 0:
    CloseSFAC()
