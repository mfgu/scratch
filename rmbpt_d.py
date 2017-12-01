#coding=gbk
from pfac.fac import *
from time import time
import sys
import math

InitializeMPI(int(sys.argv[5]))
number = int(sys.argv[1])
z = number
nele = int(sys.argv[2])
i = int(sys.argv[3])
nmax = int(sys.argv[4])

if ( nele != 2 ):
   print "核外电子数不能被此脚本处理，错误"
   exit()

#ConvertToSFAC('d%d.sf'%(i))

asym = ATOMICSYMBOL[number]

#n1 = [range(2,8),range(8,13),range(13,18),[20,23,27,31,35],[40,48,58,70],[85,100,125,150,200]]
#n2 = range(8)+[8,13,20,28,39,50,80]

n1 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20, 22, 24, 28, 32, 38, 43, 50, 65, 80, 100, 125, 150, 175, 200]
"""
if ( i == 0 ):
   n1 = [1, 2, 3, 4, 5, 6, 7]
elif ( i == 1 ):
   n1 = [8, 9, 10, 11]
elif ( i == 2 ):
   n1 = [12, 13, 15, 17, 20]
else:
   n1 = [24, 28, 32, 38, 43, 50, 65, 80, 100, 125]
"""
n2 = range(8)+[8, 13, 20, 28, 38, 50, 65, 80, 100]

n1_len = len(n1)

t0 = time()
prar=[3.057, 3.061, 3.122, 3.189, 3.261, 3.365, 3.427, 3.435, 3.476, 3.544, 3.591, 3.599, 3.642, 3.706, 3.737, 3.788, 3.775, 3.882, 3.929, 3.997, 4.074, 4.097, 4.140, 4.163, 4.188, 4.203, 4.220, 4.242, 4.270, 4.324, 4.409, 4.424, 4.482, 4.494, 4.532, 4.544, 4.614, 4.617, 4.654, 4.680, 4.743, 4.750, 4.787, 4.804, 4.839, 4.855, 4.877, 4.892, 4.912, 4.962, 5.084, 5.113, 5.162, 5.060, 5.221, 5.202, 5.251, 5.226, 5.312, 5.370, 5.342, 5.351, 5.367, 5.339, 5.413, 5.402, 5.428, 5.436, 5.463, 5.476, 5.501, 5.521, 5.526, 5.539, 5.655, 5.658, 5.684, 5.670, 5.710, 5.700, 5.851, 5.744, 5.864, 5.905, 5.815, 5.815, 5.843, 5.850, 5.857]

if ( 12 <= number <= 100):
   rp = prar[number-12]
   SetAtom(asym, -1, -1, rp)
else:
   SetAtom(asym)

eps_tmp = 1.0e-5

km1= nele - 1
km2= nele - 2
km3= nele - 3

i0 = 0
g = []

if ( nele == 2 ):
   for n in range(1, nmax+1):
       if ( n == 1 ):
          Config('g%d'%n, '%d*2'%n)
       else:
          Config('g%d'%n, '1*1 %d*1'%n)
       g.append('g%d'%n)

if (i == 0):
   p0 = 'Z%02d_ne%02d'%(number,nele)
   ListConfig(p0+'a.cfg')
if (i >= 0):
   p = 'Z%02d_ne%02di%02d'%(number,nele,i)
else:
   p = 'Z%02d_ne%02d'%(number,nele)

if (i >= 0):
   OptimizeRadial(g[i0])
   SetBoundary(nmax, eps_tmp, 1E30)
   r = GetBoundary()
   ReinitRadial(0)
   SetRadialGrid(3000, 1.1, -r[2]*1.00001, 0.0)

if (i >= 0):
   SetVP(103)
   SetMS(3, 3)
   SetSE(-1, 61)
   SetBreit(-1, 1, nmax, 1.0e-5)

SetPotentialMode(10)
OptimizeRadial(g[i0])

n0 = len(g)

if (i >= 0):
   SetBoundary(nmax, eps_tmp, 1E30)

TransitionMBPT(3, 3)
if (i >= 0):
   StructureMBPT(0,0,eps_tmp) # The 5th form of the StructureMBPT, can test 3rd order corrections
   StructureMBPT(p+'b.en', p+'b.ham', g, 
                  [-1,-1,-1,-1,-1,-1,-1,-1,15,15,15,15],
                  n1, n2, n0)  # The 1st form of the StructureMBPT
   MemENTable(p+'b.en')
   PrintTable(p+'b.en', p+'a.en')
else:
   SetTransitionCut(0.0)
   SetMixCut(0.0)
   SetAngZCut(0.0)
   h = [p+'i%02db.ham'%x for x in range(1)]
   StructureMBPT(1) # The 4th form of the StructureMBPT, can be used to indicate the levels which will be corrected with MBPT
   StructureMBPT(p+'b.en', p+'a.ham', h, g, n0) # The 3rd form of the StructureMBPT
   MemENTable(p+'b.en')
   PrintTable(p+'b.en', p+'a.lev',1)
   BasisTable(p+'a.bas')
   nlev = LevelInfor(p+'b.en', -1001)
   print 'nlev = %d'%nlev
   if ( nlev > 1000 ):
      nlev = 1000
   for m in [-1, 1, -2, 2]:
       for iU in range(1, nlev):
           if (iU/10*10 == iU ):
              print 'm=%d,  iU = %d'%(m, iU)
           for iL in range(0, iU):
               TRTable(p+'b.tr', [iL], [iU], m)
   PrintTable(p+'b.tr', p+'a.tr',1)
   """
   ztmp = int((nele - 1.0) * 5.0 / 6.0 *100)/100.0
   zeff = z - ztmp
   zeff2 = zeff*zeff*13.6057
   u = [0.002, 0.008, 0.03, 0.08, 0.20, 0.80]

   e = []
   nEi = len(u)
   for i in range(nEi):
       Ef = u[i]*zeff2
       e.append(Ef)
   #print(e)
   e0 = e

   u = [10, 30, 100]
   e = []
   nEi = len(u)
   for i in range(nEi):
       Ef = u[i]*zeff2
       e.append(Ef)
   ep = e
   eborn = -1000 * zeff2

   Lmax = 100
   Lcb = 100
   SetCEBorn(eborn)
   SetCELQR(5)
   SetCELMax(Lmax)
   SetCELCB(Lcb)
   SetCEGrid(e0)
   CETable(p+'b.ce', range(nlev), range(nlev))

   SetCEBorn(eborn,0)
   SetCEGrid(ep)
   CETable(p+'b.ce', range(nlev), range(nlev))

   PrintTable(p+'b.ce', p+'a.ce',1)
   """
FinalizeMPI()
#CloseSFAC()
