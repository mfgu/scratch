#ifndef _CFVHONE_
#define _CFVHONE_

#include "params.h"

#define NGHOST  6
#define NGHOST2 (2*NGHOST)
#define G 6.67E-8
#define MP 1.67E-24
#define PC 3.08E18
#define MS 1.9889E33
#define LS 3.9E33
#define PI 3.14159
#define C 3E10
#define KB 1.38E-16
#define SIGT 6.65E-25
#define YEAR 3.15E7
#define MMW 0.7
#define XMIN (-8.0)
#define XMAX 10.0
#define TMIN 0.0
#define TMAX 8.0

#define Min(a,b) ((a)<(b)?(a):(b))
#define Max(a,b) ((a)<(b)?(b):(a))

typedef struct _HCF_ {
  int nx, nt;
  double *xi, *te;
  double **htot, **ctot, **fion, **mwgt, *ete;
  double **fhm, **racc1, **racc2, **racc3, **racc;
  double tmin;
} HCF;

typedef struct _VH1State_ {
  double gamma, gamma1;
  double mb;
  double mdot, mdoto, mtoti, mtoto;
  double reff, reffa, reffb, reffc, reffd, refft, reffs;
  double fhm;
  double gconst;
  double ulength, utime, umass, uenergy, ulum, uvel, uden, uvol;
  double lum, lumedd, emb, medd, lumts, lumi;
  int nr, nrm, nri, nro, nrt, ilog, nleft, nright, selfg;
  double rmin, rmax, ris, ros, drg;
  double *rg, *rc, *rd, *nrc;
  double *d0, *p0, *u0, *ws;
  double di0, di, dib, vi, ti, pi, cs;
  double *d, *p, *u, *te;
  double *f0, *f1, *g0, *g1, *xi, *lumr, *lumrs, *lumc, *hr, *rf;
  double *dsrc, *gsrc, *esrc, *fion;
  double time, endtime, dt;
  long long ncycle;
  int endcycle;
  double mh, ah;
  double ri, rb, tcycle;
  HCF *hcf;
  int nmdotmax, nmdotmin;
  double vmdot, tmdotmin, tmdotmax, wmdot;
  int ncheck;
  char odir[BUFLEN];
  char fmdot[BUFLEN];
  char funit[BUFLEN];
  char fsolution[BUFLEN];
  char fhcf[BUFLEN];
  char initfile[BUFLEN];
  double tmin, hidx;
  double gsoft;
  double loc_rad;
  double lum_updm;
  int lum_updc;
  double mgas0, mgas, mgasm;
  int teq, miter;
  double dtfv, dtfc;
} VH1State;

void set_params(VH1State *vh1s);
HCF *init_hcf(char *fn, double tmin, double hidx);
void alloc_vh1state(VH1State *vh1s);
void set_units(VH1State *vh1s);
void init_vh1state(VH1State *vh1s);
double halo_gravity(double r, double m, double a);
void set_ghalo(VH1State *vh1s);
void set_gsrc0(VH1State *vh1s);
void set_gsrc1(VH1State *vh1s);
double qtemp(HCF *hc, double q, double uvel, double logx);
double interpolx(double x, int n, double *xg, int *i0, int *i1);
double interpol(double x, int n, double *xg, double *yg);
double interpol2d(double x, double y, 
		  int nx, double *xg, 
		  int ny, double *yg,
		  double **z);
double net_hc(HCF *hcf, double logxi, double logt);
double mwgt(HCF *hcf, double logxi, double logt);
double racc_line(HCF *hcf, double logxi, double logt);
double racc_abs(HCF *hcf, double logxi, double logt);
double racc_sct(HCF *hcf, double logxi, double logt);
double racc_tot(HCF *hcf, double logxi, double logt);
double ion_frac(HCF *hcf, double logxi, double logt);
void set_esrc0(VH1State *vh1s);
void set_lumr(VH1State *vh1s);
void set_esrc1(VH1State *vh1s);
void mdot_stat(VH1State *vh1s);
void evolve_vh1state(VH1State *vh1s);

double nc_scale(double m);
double mdot_scale(double n, double m);
double ri_scale(double n, double m);
double rb_scale(double m);
double tc_scale(double n, double m);
int esrc1fun(double t, const double y[], double f[], void *p);

#endif
