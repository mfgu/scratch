#coding=gbk
from pfac.fac import *
from time import time
import sys
import math

ConvertToSFAC('dn.sf')
InitializeMPI(int(sys.argv[5]))
number = int(sys.argv[1])
z = number
nele = int(sys.argv[2])
i = int(sys.argv[3])
nmax = int(sys.argv[4])

if ( nele < 3 or nele > 10 ):
   print "核外电子数不能被此脚本处理，错误"
   exit()


asym = ATOMICSYMBOL[number]

#n1 = [range(2,8),range(8,13),range(13,18),[20,23,27,31,35],[40,48,58,70],[85,100,125,150,200]]
#n2 = range(8)+[8,13,20,28,39,50,80]

n1 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20, 22, 24, 28, 32, 38, 43, 50, 65, 80, 100, 125, 150, 175, 200, 225, 250]
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
n2 = range(8)+[8, 13, 20, 28, 38, 50, 65, 80, 100, 125, 150]

n1_len = len(n1)

t0 = time()

#Z=70, r different PRA2013: 5.317; PRA2015: 5.237
#Z=80, r different PRA2013: 5.467; PRA2015: 5.475

prar2013=[3.428, 3.407, 3.476, 3.542, 3.599, 3.602, 3.612, 3.705, 3.736, 3.782, 3.776, 3.898, 3.955, 3.998, 4.079, 4.104, 4.171, 4.156, 4.23 , 4.245, 4.242, 4.244, 4.273, 4.318, 4.415, 4.41 , 4.475, 4.502, 4.526, 4.542, 4.613, 4.619, 4.655, 4.704, 4.804, 4.752, 4.787, 4.807, 4.84 , 4.855, 4.877, 4.893, 4.915, 4.962, 5.031, 5.041, 5.089, 5.099, 5.083, 5.21, 5.123, 5.192, 5.317, 5.246, 5.29, 5.299, 5.359, 5.351, 5.376, 5.401, 5.418, 5.437, 5.467, 5.483, 5.505, 5.531, 5.539, 5.578, 5.632, 5.64, 5.663, 5.67, 5.802, 5.7, 5.86]

prar2015=[3.955, 3.998, 4.079, 4.104, 4.171, 4.156, 4.230, 4.245, 4.242, 4.244, 4.273, 4.318, 4.415, 4.410, 4.475, 4.502, 4.526, 4.542, 4.613, 4.619, 4.655, 4.704, 4.804, 4.752, 4.826, 4.807, 4.840, 4.855, 4.877, 4.893, 4.915, 4.962, 5.031, 5.041, 5.089, 5.099, 5.083, 5.210, 5.123, 5.192, 5.237, 5.246, 5.290, 5.299, 5.359, 5.351, 5.376, 5.401, 5.418, 5.437, 5.475, 5.483, 5.505, 5.531, 5.539, 5.578, 5.632, 5.640, 5.663, 5.670, 5.804, 5.700, 5.861, 5.744, 5.794, 5.787, 5.816, 5.816, 5.844, 5.865, 5.886]

#if ( 18 <= number <= 92):
#   rp = prar2013[number-18]
#   SetAtom(asym, -1, -1, rp)
if ( 30 <= number <= 100):
   rp = prar2015[number-30]
   SetAtom(asym, -1, -1, rp)
else:
   SetAtom(asym)

eps_tmp = 1.0e-5

i0 = 0
g = []
km2 = nele-2
km3 = nele-3
for n in range(2, nmax+1):
    if ( n == 2 ):
       Config('g%d'%n, '1*2 %d*%d'%(n, km2))
    else:
       if ( km3 == 0 ):  
          Config('g%d'%n, '1*2 %d*1'%n)
       else:
          Config('g%d'%n, '1*2 2*%d %d*1'%(km3, n))
    g.append('g%d'%n)

if (i >= 0):
   p = 'Z%02d_ne%02di%02d'%(number,nele,i)
else:
   p = 'Z%02d_ne%02d'%(number,nele)

if (i >= 0):
   OptimizeRadial(g[i0])
   SetBoundary(nmax, eps_tmp, 1E30)
   r = GetBoundary()
   ReinitRadial(0)
   SetRadialGrid(3000, 1.1, -1e30, 0.0)
if ( nmax < 5 ):
   nbreit=5
else:
   nbreit=nmax

if (i >= 0):
   SetVP(103)
   SetMS(3, 3)
   SetSE(-1, 61)
   SetBreit(-1, 1, nbreit)

#SetPotentialMode(10)
SetPotentialMode(20,1e39,0)
OptimizeRadial(g[i0])

if (i == 0):
   p0 = 'Z%02d_ne%02d'%(number,nele)
   ListConfig(p0+'a.cfg')

n0 = len(g)

if (i >= 0):
   SetBoundary(nmax, eps_tmp, 1E30)

if (i >= 0):
   TransitionMBPT(0, 0)
   StructureMBPT(0,0,eps_tmp) # The 5th form of the StructureMBPT, can test 3rd order corrections
   StructureMBPT(p+'b.en', p+'b.ham', g, 
                  [-1,-1,-1,-1,-1,-1,-1,-1,20,20,20,20],
                  n1, n2, n0)  # The 1st form of the StructureMBPT
   MemENTable(p+'b.en')
   PrintTable(p+'b.en', p+'a.en')
else:
   TransitionMBPT(3, 3)
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
CloseSFAC()
