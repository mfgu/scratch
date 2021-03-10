from pfac.rfac import *
from pfac import fac
from pylab import *
import pandas as pd

def ltepop(ks, rs, tk, rho, eps=1e-8, miter=1024,
           ein=[5.52,10.7,22.14,40.4]):
    te = tk*8.6e-5
    z = int(rs[0].z)
    mass = fac.ATOMICMASS[z]
    nt = rho/(mass*1.67e-24)
    nk = len(ks)
    pf = zeros(nk)
    e0 = zeros(nk)
    g0 = zeros(nk)
    sf = zeros(nk)
    ei = zeros(nk)
    for ki in range(nk):
        k = ks[ki]
        r = rs[ki]
        bf = (r.j+1.0)/(r.j[0]+1.0)*exp(-(r.e-r.e0)/te)        
        e0[ki] = r.e0
        pf[ki] = sum(bf)
        g0[ki] = r.j[0]+1.0
        if ki > 0:
            ei[ki-1] = e0[ki]-e0[ki-1]
            if ki-1 < len(ein):
                ei[ki-1] = ein[ki-1]
            sf[ki] = 2*g0[ki]/g0[ki-1]*exp(-ei[ki-1]/te)
    if (nk == 1):
        return nt,nt*ks[0],array([nt/pf[0]]),pf
    
    sf *= 3.0185355e21*te**1.5
    np = zeros(nk)
    nx = nt*0.5
    for iter in range(miter):
        np[0] = 1.0
        ny = ks[0]*pf[0]*np[0]        
        for ki in range(1,nk):
            np[ki] = sf[ki]*np[ki-1]/nx
            ny += ks[ki]*pf[ki]*np[ki]
        n0 = nt/(sum(np*pf))
        np *= n0
        ny *= n0
        #print('%5d %g %g %g'%(iter, nx, ny, n0))
        if (abs(1-ny/nx) < eps):
            break
        nx = ny
    #print(np)
    return nt,nx,np,pf

def kpm(ks, rs, ts, tk, rho, md=0):
    te = tk*8.6e-5
    nt,nx,np,pf = ltepop(ks, rs, tk, rho)
    yt = 0.0
    ay = 1.83e23/fac.ATOMICMASS[int(rs[0].z)]
    at = 0.285
    for ik in range(len(ks)):
        k = ks[ik]
        r = rs[ik]
        t = ts[ik]
        i0 = int32(t[2])
        p = exp(-(r.e-r.e0)/te)*(r.j+1.0)/(r.j[0]+1.0)*np[ik]        
        xs = t[4]
        ys0 = t[5]/(t[3]+1)*p[i0]
        ys = ys0*ay*(xs**3)/(exp(xs/te)-1.0)/(tk)**4/nt
        if md == 1:
            tau = at*ys0/xs
            w = where(tau > 0.01)
            ys[w] *= (1-exp(-tau[w]))/tau[w]
        yi = sum(ys)
        yt += yi
    return yt

def kpi(x, y, tk):
    xe = 12.3984e3/x
    te = tk*8.6e-5
    a = 3.46e19/tk**4
    return sum(y*xe**3/(x*(exp(xe/te)-1)))*a*(x[1]-x[0])*2/(x[0]+x[1])

def knu(ks, rs, ts, tk, rho, lam0, lam1, dlam=0.01, tej=1.0):
    te = tk*8.6e-5
    nt,nx,np,pf = ltepop(ks, rs, tk, rho)
    xlam0 = log(lam0)
    xlam1 = log(lam1)
    xlam = arange(xlam0, xlam1+0.5*dlam, dlam)
    lam = exp(xlam)
    yr = zeros(len(lam))
    y1 = yr.copy()
    ay = 6.63e7/fac.ATOMICMASS[int(rs[0].z)]
    at = 0.285
    for ik in range(len(ks)):
        k = ks[ik]
        r = rs[ik]
        t = ts[ik]
        i0 = int32(t[2])
        p = exp(-(r.e-r.e0)/te)*(r.j+1.0)/(r.j[0]+1.0)*np[ik]   
        xs = t[4]
        de = dlam*xs
        xs = log(12.4e3/xs)
        ix = int32((xs-xlam0)/dlam)
        ys0 = t[5]/(t[3]+1)*p[i0]
        ys = ay*ys0/de/nt          
        w = where((xs >= xlam0)&(xs < xlam1))
        ys1 = ys.copy()
        tau = at*tej*ys0/t[4]
        w1 = where(tau > 0.01)
        ys1[w1] *= (1-exp(-tau[w1]))/tau[w1]
        df = pd.DataFrame({'idx': ix[w], 'data': ys[w]},
                          columns=['idx','data'])
        dfs = df.groupby('idx').sum()
        yr[dfs.index] += dfs.data
        df = pd.DataFrame({'idx': ix[w], 'data': ys1[w]},
                          columns=['idx','data'])
        dfs = df.groupby('idx').sum()
        y1[dfs.index] += dfs.data
        
    return lam,yr,y1

def fit_exop(ks, rs, ts,
             tk0=3.0, tk1=5.0, dtk=0.25,
             t0=-6.0, t1=9.0, dt=1.0,
             x0=1e3, x1=1e5, dx=0.01,
             sav=None):
    rho = 1e-13
    tk = 10**arange(tk0, tk1+0.1*dtk, dtk)
    ntk = len(tk)
    tej = 10**arange(t0, t1+0.1*dt, dt)
    nt = len(tej)
    
    a = zeros((ntk,nt,2))
    a[:,-1,1] = 1.0
    ys = None
    for k in range(ntk):
        for i in range(nt):
            x,yb,y = knu(ks, rs, ts, tk[k], rho, x0, x1, dlam=dx, tej=tej[i])
            nx = len(x)
            if ys is None:
                ys = zeros((ntk,nt,nx))
                ybs = zeros((ntk,nx))
            if i == 0:
                xd = log10(yb)                
                ybs[k] = yb
            else:
                yd = log10(ys[k,i-1]/y)/log10(tej[i]/tej[i-1])
                w = where((y > 0)&(ys[k,i-1]>0))
                a[k,i-1] = polyfit(xd[w], yd[w], 1)
                print([k,i,a[k,i-1]])
            ys[k,i] = y
        if len(ks) == 1 and not sav is None:
            fn = '%s%s%02dt%d.dat'%(sav,rs[0].asym,ks[0],k)
            with open(fn, 'w') as f:
                s = '#atom: %s\n'%rs[0].asym
                f.write(s)
                s = '#z: %d\n'%rs[0].z
                f.write(s)
                s = '#cs: %d\n'%ks[0]
                f.write(s)
                s = '#ip: %15.8E\n'%(rs[0].ei)
                f.write(s)
                s = '#g0: %d\n'%(rs[0].j[0]+1)
                f.write(s)
                s = '#e0: %15.8E\n'%rs[0].e0
                f.write(s)
                te = tk[k]*8.6e-5
                pf = sum(exp(-(rs[0].e-rs[0].e0)/te)*(rs[0].j+1)/(rs[0].j[0]+1))
                s = '#pf: %15.8E\n'%pf
                f.write(s)
                s = '#itk: %d\n'%k
                f.write(s)
                s = '#tk: %15.8E\n'%tk[k]
                f.write(s)
                s = '#rho: 1E-13\n'
                f.write(s)
                s = '#logt0: %15.8E\n'%t0
                f.write(s)
                s = '#logt1: %15.8E\n'%t1
                f.write(s)
                s = '#dlogt: %15.8E\n'%dt
                f.write(s)
                s = '#nt: %d\n'%nt
                f.write(s)
                s = '#logtk0: %15.8E\n'%tk0
                f.write(s)
                s = '#logtk1: %15.8E\n'%tk1
                f.write(s)
                s = '#dlogtk: %15.8E\n'%dtk
                f.write(s)
                s = '#ntk: %d\n'%ntk
                f.write(s)
                s = '#lam0: %15.8E\n'%x0
                f.write(s)
                s = '#lam1: %15.8E\n'%x1
                f.write(s)
                s = '#dlam: %15.8E\n'%dx                
                f.write(s)
                s = '#nlam: %d\n'%nx
                s = '#a0:'
                for m in range(nt):
                    s += ' %12.5E'%a[k,m,0]
                s += '\n'
                f.write(s)
                s = '#a1:'
                for m in range(nt):
                    s += ' %12.5E'%a[k,m,1]
                s += '\n'
                f.write(s)
                for j in range(nx):
                    s = '%12.5E %12.5E '%(x[j], ybs[k,j])
                    for m in range(nt):
                        s += ' %12.5E'%ys[k,m,j]
                    s += '\n'
                    f.write(s)
                    
    return tk,tej,x,ybs,ys,a

def exop(tk, tr, x, yb, yr, a, tk0, rho, tej):
    k0 = argmin(abs(tk-tk0))
    if (k0 == 0):
        k1 = k0+1
    elif k0 == len(tk)-1:
        k1 = k0
        k0 = k0-1
    else:
        if tk0 > tk[k0]:
            k1 = k0+1
        else:
            k1 = k0
            k0 = k0-1
    tej = tej*rho/1e-13
    w = argmin(abs(tr-tej))
    if tej < tr[w]:        
        w0 = w-1
        w1 = w
        if w0 < 0:
            w0 = w
            w1 = w0
    else:
        w0 = w
        w1 = w0+1
        if w == len(tr)-1:
            w1 = w0
    if w1 > w0:
        wf = (log(tej)-log(tr[w0]))/(log(tr[w1])-log(tr[w0]))
    else:
        wf = 0.0
    a0 = a[k0,w0]
    if w0 == len(tr)-1:
        a0 = [0.0, 1.0]
    c = zeros(len(yb[k0]))
    k = where(yb[k0]>0)
    c[k] = a0[0]*log10(yb[k0][k]) + a0[1]
    yf0 = yr[k0,w0]/(tej/tr[w0])**c
    a0 = a[k0,w1]
    if w1 == len(tr)-1:
        a0 = [0.0, 1.0]
    c = zeros(len(yb[k0]))
    c[k] = a0[0]*log10(yb[k0][k]) + a0[1]    
    yf0 = (1-wf)*yf0 + wf*yr[k0,w1]/(tej/tr[w1])**c
    
    a1 = a[k1,w0]
    if w0 == len(tr)-1:
        a1 = [0.0,1.0]
    c = zeros(len(yb[k1]))
    k = where(yb[k1]>0)
    c[k] = a1[0]*log10(yb[k1][k]) + a1[1]
    yf1 = yr[k1,w0]/(tej/tr[w0])**c
    a1 = a[k1,w1]
    if w1 == len(tr)-1:
        a1 = [0.0, 1.0]
    c = zeros(len(yb[k1]))
    c[k] = a1[0]*log10(yb[k1][k]) + a1[1]
    yf1 = (1-wf)*yf1 + wf*yr[k1,w1]/(tej/tr[w1])**c
    
    xt0 = log(tk[k0])
    xt1 = log(tk[k1])
    xt = log(tk0)
    xf = (xt-xt0)/(xt1-xt0)
    yf = zeros(len(yf0))
    w = where((yf0>0)&(yf1>0))
    yf[w] = exp((1-xf)*log(yf0[w]) + xf*log(yf1[w]))
    return yf

def mix_exop(ks, rs, a, tk, rho, tej, pkm=0):
    p = ltepop(ks, rs, tk, rho)
    p = p[2]*p[3]/p[0]
    nx = len(a[0][2])
    yf = zeros(nx)
    for k in range(len(ks)):
        yf += p[k]*exop(*a[k], tk, rho*p[k], tej)
    if pkm:
        yf = kpi(a[0][2], yf, tk)
    return yf

def cmp_exop(ks, rs, ts, tk, rho, a):
    clf()
    tr = a[1]
    for i in range(len(tr)):
        if i < len(tr)-1:
            t = sqrt(tr[i]*tr[i+1])
        else:
            t = tr[i]*(tr[i]/tr[i-1])
        x, yb, yr = knu(ks, rs, ts, tk, rho, 1e3, 1e5, tej=t)
        yf = exop(*a, tk, rho, t)
        loglog(x, yr, marker='.', linestyle='none', color='k')
        loglog(x, yf, color='r', linewidth=0.5)
        
def rdf(k0, k1):
    rs = []
    ks = []
    ts = []
    for k in range(k0, k1+1):
        r = FLEV('ta%d.en'%k)
        t = load_fac('ta%d.tr'%k)
        ks.append(k)
        rs.append(r)
        ts.append(t)
    return ks, rs, ts

