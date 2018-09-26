from numpy import *
import sys

m = int(sys.argv[1])

e0 = 300.0
nn = 1000
sig = 2.13
de = sig/2
e = [e0 + i*de for i in range(nn)]

f = open('res.txt', 'w')
for i in range(nn):
    fn = 'pol%de%03d.txt'%(m, i)
    print(fn)
    r = transpose(loadtxt(fn))
    if len(r.shape) == 1:
        r = r.reshape((8,1))
    w0 = where(logical_and(r[4] > 700, r[4] < 750))
    w1 = where(logical_and(r[4] > 800, r[4] < 850))
    y0 = 0.0
    y1 = 0.0
    p0 = 0.0
    p1 = 0.0
    s0 = 0.0
    s1 = 0.0
    if len(w0[0])>0:
        y0 = sum(r[5][w0]*r[6][w0])
        s0 = sum(r[5][w0])
    if (y0 > 0):
        p0 = sum(r[5][w0]*r[7][w0])/sum(r[5][w0])
    if len(w1[0]) > 0:
        y1 = sum(r[5][w1]*r[6][w1])
        s1 = sum(r[5][w1])
    if (y1 > 0):
        p1 = sum(r[5][w1]*r[7][w1])/sum(r[5][w1])
    f.write('%12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n'%(e[i], y0, y1, p0, p1, s0, s1))
f.close()
