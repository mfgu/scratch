from optparse import OptionParser
from pfac.fac import *
from time import time
import math
import os

t0 = time()

#parse cmdline args
p = OptionParser()
p.add_option('-z', '--z', dest='z', type='int',
             default=26, help='atomic number')
p.add_option('-n', '--n', dest='n', type='int',
             default=1, help='number of M-shell electrons')
p.add_option('-p', '--np', dest='np', type='int',
             default=1, help='number of processors')
p.add_option('-m', '--nmax', dest='nmax', type='int',
             default=5, help='max principle quantum number')
p.add_option('-i', '--imax', dest='imax', type='int',
             default=0, help='inner shell excitation')
p.add_option('-v', '--vmax', dest='vmax', type='int',
             default=-8, help='virtual orbital n max')
p.add_option('--n2max', dest='n2max', type='int',
             default=0, help='max n of the 2nd virtual orbital')
p.add_option('-r', '--nr', dest='nr', type='int',
             default=0, help='job id for the n1 split run')
p.add_option('-s', '--nsp', dest='nsp', type='int',
             default=16, help='number of n1s in split jobs')
p.add_option('-t', '--ntr', dest='ntr', type='int',
             default=0, help='mbpt transition rate option')
p.add_option('-k', '--nk', dest='nk', type='int',
             default=-1, help='k shell max excitation')
p.add_option('-c', '--csf', dest='csf', type='int',
             default=0, help='convert to sfac input file')
p.add_option('-a', '--mmax', dest='mm', type='float',
             default=0, help='maximum memory usage for radial integra caching')
p.add_option('-d', '--dry', dest='dry', type='int',
             default=0, help='dry run, stop after configuration setup')
p.add_option('--m3d', dest='m3d', type='int',
             default=0, help='max n of double excitation of M-shell')
p.add_option('--m3i', dest='m3i', type='int',
             default=0, help='max n of double excitation of L-shell')
p.add_option('--oc', dest='oc', type='string',
             default='ic', help='configuration for potential optimization')
p.add_option('--om', dest='om', type='int',
             default=20, help='optimizaiton mode')
p.add_option('--od', dest='od', type='string',
             default='', help='output directory')
p.add_option('--mcc', dest='mcc', type='int',
             default=9, help='maximum n of correlation config')
p.add_option('--mcc2', dest='mcc2', type='int',
             default=9, help='maximum n of correlation config')
p.add_option('--kmax', dest='kmax', type='int',
             default=8, help='max orbital partial wave')
p.add_option('--kcc', dest='kcc', type='int',
             default=8, help='maximum l of correlation config')
p.add_option('--acc', dest='acc', type='float',
             default=0.05, help='correlation mixing threshold')
p.add_option('--acc2', dest='acc2', type='float',
             default=0.05, help='correlation mixing threshold')
p.add_option('--hiter', dest='hiter', type='int',
             default=0, help='hamilton iteration for perturbing configs')
p.add_option('--piter', dest='piter', type='int',
             default=50, help='perturb config iteration to enlarging ci space')
p.add_option('--ptol', dest='ptol', type='float',
             default=0.01, help='perturb config cutoff tolerance')
p.add_option('--expdim', dest='expdim', type='float',
             default=0.05, help='perturb config cutoff tolerance')
p.add_option('--expdimz', dest='expdimz', type='float',
             default=1e-4, help='perturb config cutoff tolerance')
p.add_option('--mcut0', dest='mcut0', type='float',
             default=1e-3, help='correlation mixing threshold')
p.add_option('--mcut1', dest='mcut1', type='float',
             default=1e-1, help='correlation mixing threshold')
p.add_option('--mcut2', dest='mcut2', type='float',
             default=1e-1, help='correlation mixing threshold')
p.add_option('--mcut3', dest='mcut3', type='float',
             default=1.0, help='correlation mixing threshold')
p.add_option('--nrg', dest='nrg', type='int',
             default='1500', help='radial grid points')
p.add_option('--nbr', dest='nbr', type='int',
             default=-3, help='breit max n')
p.add_option('--mbr', dest='mbr', type='int',
             default=0, help='breit mode')
p.add_option('--kbr', dest='kbr', type='int',
             default=0, help='breit min n')
p.add_option('--nse', dest='nse', type='int',
             default=-2, help='self energy maxn')
p.add_option('--mse', dest='mse', type='int',
             default=41, help='self energy mode')
p.add_option('--ci', dest='ci', type='int',
             default=0, help='do ci calculation only')
p.add_option('--bas', dest='bas', type='int',
             default=1, help='print out basis and mixing coeff')
p.add_option('--pj', dest='pj', type='int',
             default=-1, help='symmetry to include')
p.add_option('--pp', dest='pp', type='int',
             default=-1, help='parity to include')
p.add_option('--jmin', dest='jmin', type='int',
             default=-1, help='min 2J to include')
p.add_option('--nj', dest='nj', type='int',
             default=1, help='number of Js to include')
p.add_option('--dnm', dest='dnm', type='int',
             default=0, help='extra delta n for config in ci')
p.add_option('--odn', dest='odn', type='int',
             default=0, help='potential boundary n')
p.add_option('--mcut', dest='mcut', type='float',
             default=0.65, help='mixing coeff cutoff for id levels')
p.add_option('--pm', dest='pm', type='int',
             default=2, help='parallel mode')
p.add_option('--pseudo', dest='pseudo', type='int',
             default=-100, help='use pseudo orbital for virtual')
p.add_option('--xdf', dest='xdf', type='float',
             default=-1, help='xdf param for pseudo orb')
p.add_option('--rand', dest='rand', type='int',
             default=11, help='randomize config list')
p.add_option('--warn', dest='warn', type='float',
             default=0.5, help='warn large mbpt terms')
p.add_option('--ignore', dest='ignore', type='float',
             default=50.0, help='ignore large mbpt terms')
p.add_option('--rc', dest='rc', type='string',
             default='', help='read config list')
p.add_option('--rh', dest='rh', type='string',
             default='', help='read hamilton')
p.add_option('--ic', dest='ic', type='string',
             default='', help='individual config')
p.add_option('--itr', dest='itr', type='int',
             default=0, help='compute CI transition rates')
p.add_option('--ice', dest='ice', type='int',
             default=0, help='compute CI excitation')

(opts, args) = p.parse_args()

print(opts)

if opts.od != '':
    x = os.system('mkdir %s'%opts.od)
odir = opts.od
if opts.ic != '':
    if odir == '':
        odir='.'
    odir = '%s/%s'%(odir,opts.ic)
    x = os.system('mkdir %s'%odir)
    
ir = opts.nr
asym = ATOMICSYMBOL[opts.z]

if opts.csf > 0:
    if ir >= 0:
        ConvertToSFAC('dl%d.sf'%ir)
    else:
        ConvertToSFAC('dl.sf')

if opts.np > 1:
    InitializeMPI(opts.np)
SetAtom(asym)

n = opts.n    
nmax = opts.nmax
pref='%s%02d'%(asym, n)
if odir != '':
    pref='%s/%s'%(odir,pref)
    
if ir >= 0 and opts.ci == 0:
    p0 = '%si%02d'%(pref,ir)
elif opts.ci > 0:
    p0 = '%sc'%pref
else:
    p0 = pref

opts.m3d = min(opts.m3d, opts.nmax)
opts.m3i = min(opts.m3i, opts.imax)

m3d=opts.m3d
if opts.rc != '':
    ReadConfig(opts.rc)

n = n-2
gc=[]
gc.append('g')
if opts.rc == '' :
    if n < 3:
        Config('g', '1s2 2s%d'%n)
    else:
        Config('g', '1s2 2s2 2p%d'%(n-2))
        
gv=['g']
v3 = []
v4 = []
for m in range(2, nmax+1):
    bc = ['2s', '2p']
    for k0 in range(2):
        for k1 in range(m):
            if m == 2 and k1 <= k0:
                continue
            gn = 'g_n2k%d_n%dk%d'%(k0,m,k1)
            gc.append(gn)
            if m == 2:
                gv.append(gn)
            elif m == 3:
                v3.append(gn)
            elif m == 4:
                v4.append(gn)
            if opts.rc == '':
                Config(1, gn, ['g'], bc[k0], m, m, k1, k1)
i=1
for m in range(2, nmax+1):
    for k0 in range(2):
        for k1 in range(m):
            if m == 2 and k1 <= k0:
                continue
            gn = 'd_n2k%d_n%dk%d'%(k0,m,k1)
            if m <= opts.m3d:
                gc.append(gn)
                #gv.append(gn)
                if opts.rc == '':
                    Config(1, gn, [gc[i]], '2*1', 2, 2)
            i = i + 1
gi=[]
i0 = len(gc)
if (opts.imax > 1):
    for m in range(2, opts.imax+1):
        bc = ['1s']
        for k0 in range(1):
            for k1 in range(m):
                gn = 'g_n1k%d_n%dk%d'%(k0, m, k1)
                gc.append(gn)
                if m == 2:
                    gi.append(gn)
                if opts.rc == '':
                    Config(1, gn, ['g'], bc[k0], m, m, k1, k1)
    i = i0
    for m in range(2, opts.imax+1):        
        for k0 in range(1):
            for k1 in range(m):
                gn = 'd_n2k%d_n%dk%d'%(k0,m,k1)
                if m <= opts.m3i:
                    gc.append(gn)
                    #gi.append(gn)
                    if opts.rc == '':
                        Config(1, gn, [gc[i]], '2*1', 2, 2)
                i = i + 1
                
if opts.vmax <= 0:
    n1 = [2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 16, 24, 36, 54, 90]
else:
    n1 = list(range(2, opts.vmax+1))
if opts.n2max <= 0:
    n2 = list(range(9))+[9, 11, 15, 23, 40]
    nn2 = len(n2)
else:
    n2 = opts.n2max
    nn2 = n2

if opts.vmax <= 0 and opts.n2max <= 0:
    nvm = n1[-1]+n2[-1]
elif opts.vmax <= 0:
    if opts.n2max < 1000:
        nvm = n1[-1]+opts.n2max
    else:
        nvm = max(n1[-1], opts.n2max/1000)
else:
    if opts.n2max < 1000:
        nvm = opts.vmax+opts.n2max
    else:
        nvm = max(opts.vmax, opts.n2max/1000)

nn = len(n1)
if opts.pm == 0:
    ncp = 0
    if (opts.nsp <= 0):
        nsp = nn
        ns = [0, nn]
    elif opts.nsp >= nn-1:
        nsp = nn-1
        ns = [0]+list(range(2,len(n1)+1))
    else:
        nsp = opts.nsp
        ns = list(range(0,nn+1,nn/nsp))
    if ns[-1] < nn:
        ns[-1] = nn
    ni = n1[ns[i]:ns[i+1]]
else:
    ncp = opts.nsp
    nsp = nn
    ns = [0, nn]
    ni = n1
    
Print('n1:', n1)
Print('nn1=%d, nn2=%d nsp=%d'%(nn, nn2, nsp))

eps = 1e-4
if (opts.nbr < -1):
    opts.nbr = max(opts.nmax, opts.imax)+abs(opts.nbr)-1
if (opts.nse < -1):
    opts.nse = max(opts.nmax, opts.imax)+abs(opts.nse)-1
#QED correction options
SetVP(103)
SetMS(3, 3)
SetSE(opts.nse, opts.mse)
SetBreit(opts.nbr, opts.mbr, -1, -1, opts.kbr)

PrintNucleus()
PrintNucleus(1, p0+'a.iso')
PrintQED()
Print('kmax=%d'%opts.kmax)
Print('ns=%d'%(len(ns)-1))
Print(ns)
if ir >= 0:
    Print(ni)
    Print(n2)
    
if opts.dry > 0:
    if ir >= 0 and opts.np > 1:
        FinalizeMPI()
    exit(0)
    
if (opts.mm > 0):
    LimitArray(-1, opts.mm)
try:    
    OptimizeRadial('g')
except:
    exit(0)
if ir >= 0 or opts.rc == '' or opts.ntr > 0:
    SetBoundary(max(nmax,opts.imax)+opts.odn, eps, 1e30)

    ReinitRadial(0)
    SetRadialGrid(opts.nrg, 1.1, -1e30, 0.0, 1.0)

    Print('opt config: %s %d'%(opts.oc, opts.om))
    SetPotentialMode(opts.om, 1e30, 0)
    try:
        if opts.oc == 'g':
            OptimizeRadial('g')
        elif opts.oc == 'gv':
            OptimizeRadial(gv)
        elif opts.oc == 'v3':
            OptimizeRadial(v3)
        elif opts.oc == 'v4':
            OptimizeRadial(v4)
        elif opts.oc == 'gv3':
            OptimizeRadial(gv+v3)
        elif opts.oc == 'gi':
            OptimizeRadial(gi)
        elif opts.oc == 'gvi':
            OptimizeRadial(gv+gi)
        elif opts.oc == 'ic' and opts.ic != '':
            OptimizeRadial(opts.ic)
    except:
        exit(0)
    
    SetBoundary(max(nmax,opts.imax)+opts.odn, eps, 1e30)
GetPotential(p0+'a.pot')

if opts.pseudo >= -2:
    Print('solve pseudo orbs')
    SolvePseudo(opts.kmax, nvm, 0, 0, opts.pseudo, opts.xdf)
    
if opts.mcut > 0:
    SetMixCut(-1, opts.mcut)

if opts.pj >= 0:
    Structure(opts.pj%2, opts.pj/2)    
elif opts.jmin >= 0 and opts.nj > 0:
    j0 = opts.jmin
    if opts.n%2  != j0%2:
        j0 += 1
    Structure(opts.pp, list(range(j0, j0+opts.nj*2, 2)))
    
Structure(opts.hiter)
Structure(-opts.piter-1, opts.ptol, opts.expdim, opts.expdimz)

ga = gc
if opts.ic != '':
    ga0 = [opts.ic]
    gc = [opts.ic]
    for c in ga:
        if c != opts.ic:
            ga0.append(c)
    ga=ga0

if opts.mcc > 1:
    ga = ga + ['cc1', 'cc2']
    if opts.rc == '':
        Config(3, 'cc1', gc, '2*1 3*1 4*1 5*1', 2, opts.mcc, 0, opts.kcc, 0, 0, opts.acc, 1)
        if opts.mcc2 > 1:
            Config(3, 'cc2', ga[len(gc):-1], '2*1 3*1 4*1 5*1', 2, opts.mcc2, 0, opts.kcc, 0, 0, -opts.acc2, gc)
    
ListConfig(p0+'a.cfg')

if opts.ci > 0:    
    Structure(p0+'b.en', p0+'b.ham', gc, ga[len(gc):], 1)
    BasisTable(p0+'a.bs')
    BasisTable(p0+'a', 10)
    MemENTable(p0+'b.en')
    PrintTable(p0+'b.en', p0+'a.en')
    TRTable(p0+'b.tr', gc[0:1], gc)
    PrintTable(p0+'b.tr', p0+'a.tr')
    exit(0)

TransitionMBPT(opts.ntr, opts.ntr)
if ir >= 0:
    StructureMBPT(opts.warn, opts.ignore)
    StructureMBPT(opts.rand, 0, opts.mcut0, opts.mcut1, opts.mcut2, opts.mcut3)
    mex = 0
    if (opts.pm == 2):
        mex = 6
    elif opts.nsp >= nn:
        mex = 1
    if opts.vmax != 0:
        StructureMBPT(abs(opts.vmax)*10+mex)
    else:
        StructureMBPT(mex)
    if opts.nk >= 0:
        StructureMBPT('1s', opts.nk)
    if opts.rh == '':
        StructureMBPT(p0+'b.en', [p0+'b.ham', p0+'b.ham0'],
                      ga, opts.kmax, ni, n2, len(gc), ncp, ir)
    else:
        StructureMBPT(p0+'b.en', [p0+'b.ham', opts.rh],
                      '', opts.kmax, ni, n2, len(gc), ncp, ir)
    if opts.bas:
        BasisTable(p0+'a.bs')
        BasisTable(p0+'a', 10)
    MemENTable(p0+'b.en')
    PrintTable(p0+'b.en', p0+'a.en')
else:
    if opts.pm == 0:
        nns = len(ns)-1
    else:
        nns = ncp
    h = [pref+'i%02db.ham'%x for x in range(nns)]
    mex = 1
    if (opts.pm == 2):
        mex += 5
    if opts.vmax != 0:
        StructureMBPT(abs(opts.vmax)*10+mex)
    else:
        StructureMBPT(mex)
    if (opts.ntr > 0):
        TransitionMBPT(p0+'b.tr', gc[0:1], gc)
    StructureMBPT(p0+'b.en', pref+'i00b.ham0', h, ga, len(gc))
    if opts.bas:
        BasisTable(p0+'a.bs')
        BasisTable(p0+'a', 10)
    MemENTable(p0+'b.en')
    PrintTable(p0+'b.en', p0+'a.en')
    if (opts.ntr > 0):
        PrintTable(p0+'b.tr', p0+'a.tr')
    if opts.itr > 0:
        TRTable(p0+'b.tr', gc, gc)
        PrintTable(p0+'b.tr', p0+'a.tr')
    if opts.ice > 0:
        CETable(p0+'b.ce', gv, gc)
        PrintTable(p0+'b.ce', p0+'a.ce')

t1 = time()
Print('all done %d %10.3E'%(opts.nr,t1-t0))
if opts.np > 1:
    FinalizeMPI()
if opts.csf > 0:
    CloseSFAC()
