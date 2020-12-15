#from pylab import *
from numpy import *
import pylab
import matplotlib.pyplot as plt
import os.path
import os
from scipy import interpolate

def lastmdot(plist, ofn):
    d = loadtxt(plist, skiprows=1, dtype=str, usecols=0, unpack=1)
    r = loadtxt(plist, skiprows=1, usecols=range(1,12), unpack=1)
    n = len(d)
    f = open(ofn, 'w')
    for i in range(n):
        u = loadtxt(d[i]+'/unit.txt', usecols=1, unpack=1)
        m = loadtxt(d[i]+'/mdot.txt', unpack=1)
        tc = 3*u[16]*3.15e7/u[0]
        w = where(m[1] >= tc)[0]
        if len(w) > 0:
            ma = (m[9][-1]*m[1][-1]-m[9][w[0]]*m[1][w[0]])/(m[1][-1]-m[1][w[0]])
        else:
            ma = m[9][-1]
        s = '%s %10.4E %10.4E %4.1f %4.1f %4.1f %4.1f %4.1f %d %10.4E'%(d[i],r[0][i],r[2][i],r[3][i],r[4][i],r[5][i],r[6][i],r[7][i],int(r[8][i]),r[9][i])
        s += ' %10.4E %10.4E %10.4E %10.4E %10.4E %10.4E'%(u[10],u[11],u[13],u[14],m[9][-1],ma)
        f.write(s+'\n')
    f.close()
    
def pmd(odir='.', t0=-1.0, t1=1e31, pt=0.001, dt=0, yr=[1e3, 1e11],
        mdmax=0, op=0, ma=0):
    u = loadtxt(odir+'/unit.txt', usecols=1, unpack=1)
    um = loadtxt(odir+'/unit.txt', usecols=0, dtype=str, unpack=1)
    w = where(um == 'tc:')[0]
    ut = u[1]/(3.15e7*u[w[0]])
    ud = 1/u[11]
    #os.system('cat mdot0.txt mdot1.txt > mdot.txt')
    d = loadtxt(odir+'/mdot.txt', unpack=1)
    
    if (t0 >= 0):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_ylim(yr[0], yr[1])
        ln = None
        w = where((d[1] >= t0/ut)&(d[1] <= t1/ut))
        ta = -2*dt
        for i in w[0]:
            f = odir+'/solution%d.txt'%(d[0][i])
            if os.path.exists(f):
                r = loadtxt(f, unpack=1)
                x = (r[3]*u[0])[6:-6]
                n = (r[5]*u[4])[6:-6]*0.75/1.67e-24
                t = r[7][6:-6]
                t1 = r[1][0]*ut
                if (t1-ta < dt):
                    continue
                ta = t1
                if ln is None:
                    ln, = ax.plot(x, n)
                    ln1, = ax.plot(x, t)
                else:
                    ln.set_ydata(n)
                    ln1.set_ydata(t)
                    #plt.draw()
                plt.pause(pt)
    else:
        if op == 0:
            pylab.clf()
        w = where(d[8]>0)
        pylab.plot(d[1][w]*ut, d[8][w]*ud)
        xt0 = d[1][1:]*ut
        yd0 = d[9][1:]*ud
        pylab.plot(xt0, yd0)
        pylab.plot(d[1][w]*ut, d[18][w])
        if ma > 0:
            fd = interpolate.interp1d(xt0, xt0*yd0, bounds_error=False,fill_value='extrapolate', kind='linear')
            nt = 2+(xt0[-1]-xt0[0])/ma
            xt = linspace(xt0[0], xt0[-1], nt)
            yt = diff(fd(xt))/diff(xt)
            pylab.plot(xt[1:], yt)
        if mdmax:
            pylab.plot(d[1][w][1:]*ut, d[11][w][1:]*ud)
        pylab.yscale('log')

