#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "cfvhone.h"
#include "params.h"

#define NWRK 2048
static double _wrk[NWRK];

PARAM params[] = {
  {PSTR, 1, "HCFile", "hcf.txt", NULL},
  {PSTR, 1, "OutDir", ".", NULL},
  {PSTR, 1, "MDotFile", "mdot.txt", NULL},
  {PSTR, 1, "UnitFile", "unit.txt", NULL},
  {PSTR, 1, "SolutionBase", "solution", NULL},
  {PINT, 1, "NMDotMin", "1", NULL},
  {PINT, 1, "NMDotMax", "50000", NULL},
  {PDBL, 1, "TMDotMin", "0.001", NULL},
  {PDBL, 1, "TMDotMax", "0.25", NULL},
  {PDBL, 1, "VMDot", "0.1", NULL},
  {PDBL, 1, "WMDot", "120", NULL},
  {PINT, 1, "NCheck", "1000", NULL},
  {PDBL, 1, "MBH", "100", NULL},
  {PDBL, 1, "MHalo", "0", NULL},
  {PDBL, 1, "RHalo", "0", NULL},
  {PDBL, 1, "RadEff", "0.3", NULL},
  {PDBL, 1, "RadEffA", "3.0", NULL},
  {PDBL, 1, "RadEffB", "10.0", NULL},
  {PDBL, 1, "RadEffC", "1.0", NULL},
  {PDBL, 1, "RadEffD", "0.1", NULL},
  {PDBL, 1, "RadEFFT", "1.0", NULL},
  {PDBL, 1, "RadEffS", "-1.0", NULL},
  {PDBL, 1, "Gamma", "1.666666666666667", NULL},
  {PDBL, 1, "LengthUnit", "0", NULL},
  {PDBL, 1, "TimeUnit", "0", NULL},
  {PDBL, 1, "MassUnit", "0", NULL},
  {PDBL, 1, "RMin", "-0.001", NULL},
  {PDBL, 1, "RMax", "-100.0", NULL},
  {PDBL, 1, "RIS", "-0.001", NULL},
  {PDBL, 1, "ROS", "-100.0", NULL},
  {PINT, 1, "NGrid", "0", NULL},
  {PDBL, 1, "DRGrid", "1.05", NULL},
  {PINT, 1, "LogScale", "1", NULL},
  {PINT, 1, "SelfGravity", "0", NULL},
  {PINT, 1, "NLeft", "1", NULL},
  {PINT, 1, "NRight", "2", NULL},
  {PDBL, 1, "EndTime", "0", NULL},
  {PINT, 1, "EndCycle", "0", NULL},
  {PDBL, 1, "InitDens", "1e5", NULL},
  {PDBL, 1, "InitDensB", "1e5", NULL},
  {PDBL, 1, "InitVel", "0.0", NULL},
  {PDBL, 1, "InitTemp", "1E4", NULL},
  {PSTR, 1, "InitFile", "0", NULL},
  {PDBL, 1, "TMin", "1E4", NULL},
  {PDBL, 1, "HeatIndex", "5", NULL},
  {PDBL, 1, "MGasMax", "0", NULL},
  {PINT, 1, "TEQ", "2", NULL},
  {PINT, 1, "MIter", "10000", NULL},
  {PDBL, 1, "GSoft", "0.0", NULL},
  {PDBL, 1, "LocRad", "0.0", NULL},
  {PDBL, 1, "LumUpdM", "1E-2", NULL},
  {PINT, 1, "LumUpdC", "250", NULL},
  {PDBL, 1, "DtFv", "1.0", NULL},
  {PDBL, 1, "DtFc", "0.75", NULL},
  {PDBL, 1, "LumTS", "1.0", NULL},
  {PNUL, 1, "", "", NULL}
};

void set_params(VH1State *vh1s) {
  int i;
   
  i = 0;
  params[i++].addr = vh1s->fhcf;
  params[i++].addr = vh1s->odir;
  params[i++].addr = vh1s->fmdot;
  params[i++].addr = vh1s->funit;
  params[i++].addr = vh1s->fsolution;
  params[i++].addr = &vh1s->nmdotmin;
  params[i++].addr = &vh1s->nmdotmax;
  params[i++].addr = &vh1s->tmdotmin;
  params[i++].addr = &vh1s->tmdotmax;
  params[i++].addr = &vh1s->vmdot;
  params[i++].addr = &vh1s->wmdot;
  params[i++].addr = &vh1s->ncheck;
  params[i++].addr = &vh1s->mb;
  params[i++].addr = &vh1s->mh;
  params[i++].addr = &vh1s->ah;
  params[i++].addr = &vh1s->reff;
  params[i++].addr = &vh1s->reffa;
  params[i++].addr = &vh1s->reffb;
  params[i++].addr = &vh1s->reffc;
  params[i++].addr = &vh1s->reffd;
  params[i++].addr = &vh1s->refft;
  params[i++].addr = &vh1s->reffs;
  params[i++].addr = &vh1s->gamma;
  params[i++].addr = &vh1s->ulength;
  params[i++].addr = &vh1s->utime;
  params[i++].addr = &vh1s->umass;
  params[i++].addr = &vh1s->rmin;
  params[i++].addr = &vh1s->rmax;
  params[i++].addr = &vh1s->ris;
  params[i++].addr = &vh1s->ros;
  params[i++].addr = &vh1s->nr;
  params[i++].addr = &vh1s->drg;
  params[i++].addr = &vh1s->ilog;
  params[i++].addr = &vh1s->selfg;
  params[i++].addr = &vh1s->nleft;
  params[i++].addr = &vh1s->nright;
  params[i++].addr = &vh1s->endtime;
  params[i++].addr = &vh1s->endcycle;
  params[i++].addr = &vh1s->di;
  params[i++].addr = &vh1s->dib;
  params[i++].addr = &vh1s->vi;
  params[i++].addr = &vh1s->ti;
  params[i++].addr = vh1s->initfile;
  params[i++].addr = &vh1s->tmin;
  params[i++].addr = &vh1s->hidx;
  params[i++].addr = &vh1s->mgasm;
  params[i++].addr = &vh1s->teq;
  params[i++].addr = &vh1s->miter;
  params[i++].addr = &vh1s->gsoft;
  params[i++].addr = &vh1s->loc_rad;
  params[i++].addr = &vh1s->lum_updm;
  params[i++].addr = &vh1s->lum_updc;
  params[i++].addr = &vh1s->dtfv;
  params[i++].addr = &vh1s->dtfc;
  params[i++].addr = &vh1s->lumts;
}

double equil_temp(int nt, double *te, double *htot, double *ctot) {
  int i;
  for (i = 0; i < nt; i++) {
    if (htot[i] < ctot[i]) break;
  }
  if (i == 0) return te[0];
  if (i == nt) return te[nt-1];
  double y0 = log(htot[i-1]/ctot[i-1]);
  double y1 = log(htot[i]/ctot[i]);
  double x0 = te[i-1];
  double x1 = te[i];
  return x0 - y0*(x1-x0)/(y1-y0);
}

HCF *init_hcf(char *fn, double tmin, double hidx) {
  HCF *hc;
  FILE *f;
  char buf[BUFLEN], *p;
  int n, m, i, j;
  double xi, te, nh, h1, h2, ne, np, htot, ctot;
  double racc1, racc2, racc3, rho;

  f = fopen(fn, "r");
  if (f == NULL) {
    printf("cannot open file: %s\n", fn);
    return NULL;
  }
  hc = malloc(sizeof(HCF));
  while (1) {
    p = fgets(buf, BUFLEN, f);    
    if (p == NULL) break;
    n = sscanf(buf, "# %d %d", &hc->nx, &hc->nt);
    if (n == 2) {
      hc->xi = malloc(sizeof(double)*hc->nx);
      hc->te = malloc(sizeof(double)*hc->nt);
      hc->htot = malloc(sizeof(double *)*hc->nx);
      hc->ctot = malloc(sizeof(double *)*hc->nx);
      hc->fion = malloc(sizeof(double *)*hc->nx);
      hc->fhm = malloc(sizeof(double *)*hc->nx);      
      hc->mwgt = malloc(sizeof(double *)*hc->nx);   
      hc->racc1 = malloc(sizeof(double *)*hc->nx);
      hc->racc2 = malloc(sizeof(double *)*hc->nx);
      hc->racc3 = malloc(sizeof(double *)*hc->nx);
      hc->racc = malloc(sizeof(double *)*hc->nx);
      hc->ete = malloc(sizeof(double)*hc->nx);
      for (i = 0; i < hc->nx; i++) {
	hc->htot[i] = malloc(sizeof(double)*hc->nt);
	hc->ctot[i] = malloc(sizeof(double)*hc->nt);
	hc->fion[i] = malloc(sizeof(double)*hc->nt);
	hc->fhm[i] = malloc(sizeof(double)*hc->nt);
	hc->mwgt[i] = malloc(sizeof(double)*hc->nt);
	hc->racc1[i] = malloc(sizeof(double)*hc->nt);
	hc->racc2[i] = malloc(sizeof(double)*hc->nt);
	hc->racc3[i] = malloc(sizeof(double)*hc->nt);
	hc->racc[i] = malloc(sizeof(double)*hc->nt);
      }
      break;
    }
  }
  m = 0;
  while (1) {
    if (NULL == fgets(buf, BUFLEN, f)) break;
    n = sscanf(buf, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
	       &xi, &te, &nh, &rho, &h1, &h2, &ne, &np, &htot, &ctot, 
	       &racc1, &racc2, &racc3);
    if (n != 13) break;
    i = m/hc->nt;
    j = m%hc->nt;
    hc->xi[i] = xi;
    hc->te[j] = te;
    hc->htot[i][j] = log10(htot/(nh*nh));
    hc->ctot[i][j] = log10(ctot/(nh*nh));
    hc->racc1[i][j] = log10(racc1/nh);
    hc->racc2[i][j] = log10(racc2/nh);
    hc->racc3[i][j] = log10(racc3/nh);
    hc->racc[i][j] = log10((racc1+racc2)/nh);
    hc->mwgt[i][j] = rho/(MP*np);
    if (h2 < 0.5) {
      hc->fion[i][j] = h2;
    } else {
      hc->fion[i][j] = 1-h1;
    }
    hc->fhm[i][j] = nh*MP/rho;
    m++;
  }
  fclose(f);
  if (hc->nt > NWRK) {
    printf("hcf tgrid number %d > NWRK %d, increase NWRK and recompile\n",
	   hc->nt, NWRK);
    exit(1);
  }
  if (hidx > 0) {
    te = log10(tmin);
    h1 = net_hc(hc, XMIN, te);
    if (h1 < 0) {
      h1 = log10(-h1);
      for (j = 0; j < hc->nt; j++) {
	h2 = pow(10, h1-hidx*(hc->te[j]-te));
	for (i = 0; i < hc->nx; i++) {
	  hc->htot[i][j] = log10(pow(10, hc->htot[i][j]) + h2);
	}
      }
    }
  }
  for (i = 0; i < hc->nx; i++) {
    hc->ete[i] = equil_temp(hc->nt, hc->te, hc->htot[i], hc->ctot[i]);
  }    
  return hc;
}

void full_fname(char *odir, char *fn) {
  char buf[BUFLEN];

  sprintf(buf, "%s/%s", odir, fn);
  strcpy(fn, buf);
}

void alloc_vh1state(VH1State *vh1s) {
  int i;

  if (strlen(vh1s->odir) > 0) {
    full_fname(vh1s->odir, vh1s->fmdot);
    full_fname(vh1s->odir, vh1s->funit);
    full_fname(vh1s->odir, vh1s->fsolution);
    full_fname(vh1s->odir, vh1s->initfile);
  }
  vh1s->gamma1 = vh1s->gamma-1;
  vh1s->di0 = vh1s->di;
  vh1s->tcycle = tc_scale(vh1s->di0, vh1s->mb);
  vh1s->rb = rb_scale(vh1s->mb);
  if (vh1s->reff > 0 ||
      (vh1s->reffa > 0 && vh1s->reffb > 0) ||
      vh1s->reffs > 0) {
    vh1s->ri = ri_scale(vh1s->di0, vh1s->mb);
  } else {
    vh1s->ri = vh1s->rb;
  }
  double r0 = Min(vh1s->ri, vh1s->rb);
  double r1 = Max(vh1s->ri, vh1s->rb);
  double dr;
  if (vh1s->rmin < 0) {
    vh1s->rmin = r0*fabs(vh1s->rmin);
  }
  if (vh1s->ris < 0) {
    vh1s->ris = r0*fabs(vh1s->ris);
  }
  if (vh1s->rmax < 0) {
    vh1s->rmax = r1*fabs(vh1s->rmax);
  } 
  if (vh1s->ros < 0) {
    vh1s->ros = r1*fabs(vh1s->ros);
  }
  if (vh1s->ilog == 1 && vh1s->drg <= 0 &&
      (vh1s->ris < vh1s->rmin || vh1s->ros > vh1s->rmax)) {
    vh1s->drg = exp((log(vh1s->rmax)-log(vh1s->rmin))/(vh1s->nr-1));
  }
  if (vh1s->drg > 0) {
    if (vh1s->ilog == 0) {
      vh1s->nrm = 1+(vh1s->rmax-vh1s->rmin)/vh1s->drg;
      vh1s->nri = 0;
      vh1s->nro = 0;
      vh1s->ris = vh1s->rmin;
      vh1s->ros = vh1s->rmax;
    } else if (vh1s->ilog == 1) {
      vh1s->nrm = 1+(log(vh1s->rmax)-log(vh1s->rmin))/log(vh1s->drg);
      if (vh1s->ris < vh1s->rmin) {
	dr = vh1s->rmin*(vh1s->drg-1.0);
	vh1s->nri = 1+(vh1s->rmin - vh1s->ris)/dr;
      } else {
	vh1s->nri = 0;
      }
      if (vh1s->ros > vh1s->rmax) {
	dr = vh1s->rmax*(vh1s->drg-1.0);
	vh1s->nro = 1+(vh1s->ros-vh1s->rmax)/dr;
      } else {
	vh1s->nro = 0;
      }
    }
  }    
  vh1s->nr = vh1s->nrm + vh1s->nri + vh1s->nro;
  vh1s->nrt = vh1s->nr + NGHOST2;
  vh1s->rg = (double *) malloc(sizeof(double)*vh1s->nrt);
  vh1s->rc = (double *) malloc(sizeof(double)*vh1s->nrt);
  vh1s->nrc = (double *) malloc(sizeof(double)*vh1s->nrt);
  vh1s->rd = (double *) malloc(sizeof(double)*vh1s->nrt);
  vh1s->d0 = (double *) malloc(sizeof(double)*vh1s->nrt);
  vh1s->p0 = (double *) malloc(sizeof(double)*vh1s->nrt);
  vh1s->u0 = (double *) malloc(sizeof(double)*vh1s->nrt);
  vh1s->d = (double *) malloc(sizeof(double)*vh1s->nrt);
  vh1s->p = (double *) malloc(sizeof(double)*vh1s->nrt);
  vh1s->u = (double *) malloc(sizeof(double)*vh1s->nrt);
  vh1s->te = (double *) malloc(sizeof(double)*vh1s->nrt);
  vh1s->gsrc = (double *) malloc(sizeof(double)*vh1s->nrt);
  vh1s->esrc = (double *) malloc(sizeof(double)*vh1s->nrt);
  vh1s->dsrc = (double *) malloc(sizeof(double)*vh1s->nrt);
  vh1s->f0 = (double *) malloc(sizeof(double)*vh1s->nrt);
  vh1s->f1 = (double *) malloc(sizeof(double)*vh1s->nrt);
  vh1s->g0 = (double *) malloc(sizeof(double)*vh1s->nrt);
  vh1s->g1 = (double *) malloc(sizeof(double)*vh1s->nrt);
  vh1s->rf = (double *) malloc(sizeof(double)*vh1s->nrt);
  vh1s->xi = (double *) malloc(sizeof(double)*vh1s->nrt);
  vh1s->lumr = (double *) malloc(sizeof(double)*vh1s->nrt);
  vh1s->lumrs = (double *) malloc(sizeof(double)*vh1s->nrt);
  vh1s->lumc = (double *) malloc(sizeof(double)*vh1s->nrt);
  vh1s->hr = (double *) malloc(sizeof(double)*vh1s->nrt);
  vh1s->fion = (double *) malloc(sizeof(double)*vh1s->nrt);
  vh1s->ws = (double *) malloc(sizeof(double)*vh1s->nrt);
  for (i = 0; i < vh1s->nrt; i++) {
    vh1s->dsrc[i] = 0.0;
    vh1s->gsrc[i] = 0.0;
    vh1s->esrc[i] = 0.0;
    vh1s->lumr[i] = 0.0;
    vh1s->lumrs[i] = 0.0;
    vh1s->lumc[i] = 0.0;
    vh1s->hr[i] = 0.0;
    vh1s->xi[i] = 0.0;
    vh1s->fion[i] = 0.0;
  }
  vh1s->mdot = 0.0;
  if (vh1s->hcf == NULL) {
    vh1s->hcf = init_hcf(vh1s->fhcf, vh1s->tmin, vh1s->hidx);
  }

  vh1s->hcf->tmin = log10(vh1s->tmin);
  vh1s->fhm = vh1s->hcf->fhm[0][0];
  vh1s->di *= MP/vh1s->fhm;
  vh1s->dib *= MP/vh1s->fhm;
  double mmw = mwgt(vh1s->hcf, XMIN, log10(vh1s->ti));
  vh1s->cs = sqrt(vh1s->gamma*KB*vh1s->ti/(mmw*MP));
  vh1s->pi = vh1s->cs*vh1s->cs*vh1s->di/vh1s->gamma;
  if (1+vh1s->ulength == 1) {
    vh1s->ulength = vh1s->rb*PC;
  } 
  if (1+vh1s->utime == 1) {
    vh1s->utime = vh1s->rb*PC/vh1s->cs;
  }
  if (1+vh1s->umass == 1) {
    vh1s->umass = vh1s->mb*MS;
  }
  set_units(vh1s);
  if (vh1s->refft > 0) {
    vh1s->refft *= vh1s->tcycle*YEAR/vh1s->utime;
  }
  if (vh1s->lumts > 0) {
    vh1s->lumts *= (vh1s->tcycle*YEAR/vh1s->utime)*(vh1s->cs/C);
  }
  vh1s->lum = 0.0;
  grid_(&vh1s->nrm, &vh1s->rmin, &vh1s->rmax,
	&vh1s->nri, &vh1s->ris,
	&vh1s->nro, &vh1s->ros,
	vh1s->rg, vh1s->rc, vh1s->rd, &vh1s->ilog);
  set_ghalo(vh1s);
  init_vh1state(vh1s);
}

void free_vh1state(VH1State *vh1s) {
  free(vh1s->rg);
  free(vh1s->rc);
  free(vh1s->nrc);
  free(vh1s->rd);
  free(vh1s->d0);
  free(vh1s->p0);
  free(vh1s->u0);
  free(vh1s->d);
  free(vh1s->p);
  free(vh1s->u);
  free(vh1s->te);
  free(vh1s->gsrc);
  free(vh1s->esrc);
  free(vh1s->dsrc);
  free(vh1s->f0);
  free(vh1s->f1);
  free(vh1s->g0);
  free(vh1s->g1);
  free(vh1s->rf);
  free(vh1s->xi);
  free(vh1s->lumr);
  free(vh1s->lumrs);
  free(vh1s->lumc);
  free(vh1s->hr);
  free(vh1s->fion);
  free(vh1s->ws);
}

void set_units(VH1State *vh1s) {
  FILE *f;

  f = fopen(vh1s->funit, "w");
  if (f == NULL) {
    printf("cannot open file: %s\n", vh1s->funit);
    exit(1);
  }
  vh1s->uvel = vh1s->ulength/vh1s->utime;
  vh1s->uenergy = vh1s->umass*vh1s->uvel*vh1s->uvel;
  vh1s->ulum = vh1s->uenergy/vh1s->utime;
  vh1s->uvol = vh1s->ulength*vh1s->ulength*vh1s->ulength;
  vh1s->uden = vh1s->umass/vh1s->uvol;
  double ug = vh1s->uvol/(vh1s->umass*vh1s->utime*vh1s->utime);
  vh1s->gconst = G/ug;

  vh1s->ris *= PC/vh1s->ulength;
  vh1s->rmin *= PC/vh1s->ulength;
  vh1s->rmax *= PC/vh1s->ulength;
  vh1s->ros *= PC/vh1s->ulength;
  vh1s->mb *= MS/vh1s->umass;
  vh1s->mh *= MS/vh1s->umass;
  if (vh1s->ah+1 == 1) {
    vh1s->ah = (2e5*PC)*pow(vh1s->mh*vh1s->umass/(1.6e14*MS), 1.0/3.0);
    vh1s->ah /= vh1s->ulength;
  } else {
    vh1s->ah *= PC/vh1s->ulength;
  }
  if (1+vh1s->endtime == 1) {
    vh1s->endtime = 1e30*vh1s->rmax/(vh1s->cs/vh1s->uvel);
  } else if (vh1s->endtime > 0) {
    vh1s->endtime *= YEAR/vh1s->utime;
  } else {
    vh1s->endtime = -vh1s->endtime*vh1s->rmax/(vh1s->cs/vh1s->uvel);
  }
  vh1s->mgasm = vh1s->mgasm*MS/vh1s->umass;
  
  double c2 = pow(C/vh1s->uvel, 2);
  vh1s->lumedd = 4.0*PI*G*(vh1s->mb*vh1s->umass)*MP*C/SIGT/vh1s->ulum;
  vh1s->medd = vh1s->lumedd/c2;

  double umd = vh1s->umass/vh1s->utime;
  double umds = umd/(MS/YEAR);
  double rb0, rb1;
  rb0 = PI*pow(G*vh1s->mb*vh1s->umass,2)*vh1s->dib/pow(vh1s->cs,3);
  rb1 = rb0*4.48*pow(vh1s->gamma,1.5);
  rb0 /= umd;
  rb1 /= umd;

  fprintf(f, "ulength: %15.8E %15.8E\n", vh1s->ulength, vh1s->ulength/PC);
  fprintf(f, "utime: %15.8E %15.8E\n", vh1s->utime, vh1s->utime/YEAR);
  fprintf(f, "umass: %15.8E %15.8E\n", vh1s->umass, vh1s->umass/MS);
  fprintf(f, "uvel: %15.8E %15.8E\n", vh1s->uvel, vh1s->uvel/1e5);
  fprintf(f, "uden: %15.8E %15.8E\n", vh1s->uden, 0.0);
  fprintf(f, "uenergy: %15.8E %15.8E\n", vh1s->uenergy, 0.0);
  fprintf(f, "ulum: %15.8E %15.8E\n", vh1s->ulum, vh1s->ulum/LS);
  fprintf(f, "uvol: %15.8E %15.8E\n", vh1s->uvol, 0.0);
  fprintf(f, "gconst: %15.8E %15.8E\n", vh1s->gconst, 0.0);
  fprintf(f, "umdot: %15.8E %15.8E\n", umd, umds);
  fprintf(f, "lumedd: %15.8E %15.8E\n", vh1s->lumedd,
	  vh1s->lumedd*vh1s->ulum/LS);
  fprintf(f, "medd: %15.8E %15.8E\n", vh1s->medd, vh1s->medd*umds);
  fprintf(f, "fhm: %15.8E %15.8E\n", vh1s->fhm, 0.0);
  fprintf(f, "abondi: %15.8E %15.8E\n", rb0, rb0*umds);
  fprintf(f, "ibondi: %15.8E %15.8E\n", rb1, rb1*umds);
  fprintf(f, "ri: %15.8E\n", vh1s->ri);
  fprintf(f, "rb: %15.8E\n", vh1s->rb);
  fprintf(f, "tc: %15.8E\n", vh1s->tcycle);
  fprintf(f, "ngrid: %d\n", vh1s->nr);
  fprintf(f, "nrm: %d\n", vh1s->nrm);
  fprintf(f, "nri: %d\n", vh1s->nri);
  fprintf(f, "nro: %d\n", vh1s->nro);
  fprintf(f, "ris: %15.8E\n", vh1s->ris);
  fprintf(f, "rmin: %15.8E\n", vh1s->rmin);
  fprintf(f, "rmax: %15.8E\n", vh1s->rmax);
  fprintf(f, "ros: %15.8E\n", vh1s->ros);
  fprintf(f, "mb: %15.8E\n", vh1s->mb);
  fprintf(f, "mh: %15.8E\n", vh1s->mh);
  fprintf(f, "ah: %15.8E\n", vh1s->ah);
  fprintf(f, "cs: %15.8E\n", vh1s->cs);
  fprintf(f, "mgasm: %15.8E\n", vh1s->mgasm);
  fprintf(f, "tend: %15.8E\n", vh1s->endtime);
  int i;
  for (i = 0; i < vh1s->hcf->nx; i++) {
    fprintf(f, "ete: %15.8E %15.8E\n", vh1s->hcf->xi[i], vh1s->hcf->ete[i]);
  }
  fflush(f);
  fclose(f);
}

void init_vh1state(VH1State *vh1s) {
  int i, j, n;
  double d, a, b;
  FILE *f;
  char buf[BUFLEN];

  f = fopen(vh1s->initfile, "r");
  if (f == NULL) {
    double cs = vh1s->cs/vh1s->uvel;
    double a2 = cs*cs;
    double ra = 2*vh1s->gconst*vh1s->mb/a2;
    a = (log(vh1s->di)-log(vh1s->dib))/(log(vh1s->rc[NGHOST])-log(vh1s->rc[NGHOST+vh1s->nr-1]));
    b = log(vh1s->di);    
    for (i = 0; i < vh1s->nrt; i++) {
      d = exp(b + a*(log(vh1s->rc[i])-log(vh1s->rc[NGHOST])));
      double x = vh1s->rc[i]/ra;
      d *= 1+0.5/x;      
      vh1s->d0[i] = d/vh1s->uden;
      vh1s->p0[i] = vh1s->pi/(vh1s->uvel*vh1s->uvel*vh1s->uden)*(d/vh1s->dib);
      if (vh1s->vi > 0.9e30) {
	vh1s->u0[i] = -0.25*cs*pow(x,-2)/(1+0.5/x);
      } else {
	vh1s->u0[i] = vh1s->vi/vh1s->uvel;
      }
      vh1s->d[i] = vh1s->d0[i];
      vh1s->p[i] = vh1s->p0[i];
      vh1s->u[i] = vh1s->u0[i];
      vh1s->te[i] = vh1s->ti;
    }
    vh1s->mgas0 = 0;
    vh1s->mgas = 0;
    for (i = NGHOST; i < vh1s->nr+NGHOST; i++) {
      vh1s->mgas0 += 4*PI*vh1s->rc[i]*vh1s->rc[i]*vh1s->rd[i]*vh1s->d[i];
    }
    vh1s->mgas = vh1s->mgas;

    for (i = 0; i < NGHOST; i++) {
      vh1s->d0[i] = vh1s->d[i] = vh1s->d0[NGHOST];
      vh1s->d0[i+vh1s->nr+NGHOST-1] = vh1s->d[i+vh1s->nr+NGHOST-1] = vh1s->d0[vh1s->nr+NGHOST-1];
      vh1s->p0[i] = vh1s->p[i] = vh1s->p0[NGHOST];
      vh1s->p0[i+vh1s->nr+NGHOST-1] = vh1s->p[i+vh1s->nr+NGHOST-1] = vh1s->p0[vh1s->nr+NGHOST-1];
    }
    vh1s->mdot = -4*PI*vh1s->rc[NGHOST]*vh1s->rc[NGHOST];
    vh1s->mdot *= vh1s->d0[NGHOST]*vh1s->u0[NGHOST];
    vh1s->time = 0.0;
    vh1s->dt = 0.0;
    vh1s->ncycle = 0;
  } else {
    //printf("initialize from solution file\n");
    i = 0;
    while (1) {
      if (NULL == fgets(buf, BUFLEN, f)) break;
      if (buf[0] == '#') {
	n = sscanf(buf, "# %ll %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		   &vh1s->ncycle, &vh1s->time, &vh1s->dt, &vh1s->mdot,
		   &vh1s->mtoti, &vh1s->mtoto,
		   &vh1s->mb, &vh1s->emb, &vh1s->lumedd,
		   &vh1s->lumi, &vh1s->lum, &vh1s->reff);
	continue;
      }
      n = sscanf(buf,
		 "%ll %lf %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
		 &vh1s->ncycle, &vh1s->time, &j, vh1s->rc+i, vh1s->p+i,
		 vh1s->d+i, vh1s->u+i, vh1s->te+i, vh1s->f0+i, vh1s->g0+i, 
		 vh1s->rf+i, vh1s->lumr+i, vh1s->xi+i,
		 vh1s->fion+i, vh1s->hr+i, &a, vh1s->lumc+i);
      if (n != 17) continue;
      if (i != j) {
	printf("init file not consistent with params\n");
	exit(1);
      }
      vh1s->d0[i] = vh1s->d[i];
      vh1s->p0[i] = vh1s->p[i];
      vh1s->u0[i] = vh1s->u[i];
      i++;
    }
  }

  init_(&vh1s->gamma, &vh1s->mb, &vh1s->nrt, vh1s->rg, vh1s->rc, vh1s->rd,
	&vh1s->ilog, vh1s->d0, vh1s->p0, vh1s->u0, &vh1s->nleft, &vh1s->nright);
}

double halo_gravity(double r, double m, double a) {
  int i;
  double x, mr, r2;

  if (1 + m == 1 || 1 + a == 1) return 0.0;
  x = r/a;        
  mr = m*(1 - (1/(1+x)) - (x/((1+x)*(1+x))));
  r2 = r*r;
  return -mr/r2;
}

void set_ghalo(VH1State *vh1s) {
  double mh;
  
  mh = vh1s->gconst * vh1s->mh;
  for (int i = 0; i < vh1s->nrt; i++) {
    vh1s->g0[i] = halo_gravity(vh1s->g0[i], mh, vh1s->ah);
  }
}

void set_gsrc0(VH1State *vh1s) {
  double mr, hdt, mr0;

  mr = vh1s->gconst * vh1s->emb;
  mr0 = vh1s->gconst * vh1s->mb/vh1s->lumedd;
  vh1s->mgas = 0;
  double gsr = 0.0;
  if (vh1s->gsoft > 0) {
    gsr = vh1s->gsoft*vh1s->ris;
    gsr *= gsr;
  }
  int i;
  double a, r0, r1, r2, mhr;
  if (vh1s->selfg) {
    for (i = 0; i < vh1s->nrt; i++) {
      if (i == 0) r0 = vh1s->rg[0];
      else r0 = vh1s->rc[i-1];
      r1 = vh1s->rc[i];
      r2 = r1*r1;
      a = (4.0*PI/3.0) * vh1s->d[i] * (r2*r1 - r0*r0*r0);
      vh1s->ws[i] = a * vh1s->gconst;       
    }
  } else {
    for (i = 0; i < vh1s->nrt; i++) {
      vh1s->ws[i] = mr;
    }
  }
  if (vh1s->selfg) {
    vh1s->ws[0] += mr;
    for (i = 1; i < vh1s->nrt; i++) {
      vh1s->ws[i] += vh1s->ws[i-1];
    }
  }  
  for (i = 0; i < vh1s->nrt; i++) {
    r1 = vh1s->rc[i];
    r2 = r1*r1+gsr;
    vh1s->f0[i] = -vh1s->ws[i]/r2;
    if (vh1s->loc_rad > 0 && vh1s->lumc[i] > 0) {
      a = mr0*vh1s->loc_rad*vh1s->lumc[i]/r2;
      vh1s->f0[i] += a;
    }
    if (vh1s->fion[i] < 1) {
      double alum = vh1s->lumr[0]-vh1s->lumr[i];
      if (alum > 0) {
	vh1s->f0[i] -= mr0*alum*(1-vh1s->fion[i])/r2;
      }
    }    
    if (vh1s->lumr[i] > 0) {
      vh1s->rf[i] = vh1s->hr[i]/(vh1s->d[i]*C/vh1s->uvel);
    } else {
      vh1s->rf[i] = 0.0;
    }
    vh1s->gsrc[i] = (vh1s->f0[i] + vh1s->g0[i] + vh1s->rf[i]);
  }
}

void set_dt(VH1State *vh1s) {
  int i;
  double ridt = 0.0;
  double xv, sv, svel = 0.0;
  for (i = NGHOST-3; i < vh1s->nrt-NGHOST+3; i++) {
    sv = sqrt(vh1s->gamma*vh1s->p[i]/vh1s->d[i])/vh1s->rd[i];
    double u = fabs(vh1s->u[i]);
    double g = fabs(vh1s->gsrc[i]);
    double gx = 2*g*vh1s->rd[i];
    double u2 = u*u;
    if (gx < 1e-5*u2) {
      xv = u/vh1s->rd[i];
    } else {
      xv = g/(sqrt(u2+gx)-u);
    }
    if (ridt < xv) ridt = xv;
    if (svel < sv) svel = sv;
  }
  ridt *= vh1s->dtfv;
  if (ridt < svel) ridt = svel;
  double dt = vh1s->dtfc/ridt;
  if (vh1s->dt > 0) {
    double dt1 = 1.2*vh1s->dt;
    if (dt > dt1) dt = dt1;
  }
  setdt_(&dt, &svel);
  vh1s->dt = dt;
  vh1s->time += dt;
  vh1s->ncycle++;
}

void set_gsrc1(VH1State *vh1s) {
  double mhr;
  
  mhr = vh1s->gconst *vh1s->mh;
  int i;
  double r1, rn;
  for (i = 0; i < vh1s->nrt; i++) {
    r1 = vh1s->rc[i];
    rn = vh1s->nrc[i];
    vh1s->f1[i] = vh1s->f0[i]*(r1*r1)/(rn*rn);
    vh1s->g1[i] = halo_gravity(rn, mhr, vh1s->ah);
    vh1s->gsrc[i] = (vh1s->f0[i] + vh1s->f1[i] + 
		     vh1s->g0[i] + vh1s->g1[i] + 2*vh1s->rf[i]);
  }
}

double qtemp(HCF *hc, double pd, double uvel, double logx) {
  double mmw, logte, te;
  int i;
  for (i = 0; i < hc->nt; i++) {
    _wrk[i] = pow(10,hc->te[i])/mwgt(hc, logx, hc->te[i]);
  }
  for (i = 0; i < hc->nt-1; i++) {
    if (_wrk[i+1]>_wrk[i]) break;
  }
  mmw = pd*uvel*uvel*MP/KB;
  logte = interpol(mmw, hc->nt-i, &(_wrk[i]), &hc->te[i]);
  te = pow(10, logte);
  return te;
}

double interpolx(double x, int n, double *xg, int *i0, int *i1) {
  double f;
  int i;
  
  if (xg[0] < xg[n-1]) {
    if (x <= xg[0]) {
      *i0 = 0;
      *i1 = 1;
    } else if (x >= xg[n-1]) {
      *i0 = n-1;
      *i1 = n-2;
    } else {
      *i0 = 0;
      *i1 = n-1;
      while (*i1 - *i0 > 1) {
	i = (*i0 + *i1)/2;
	if (x <= xg[i]) {
	  *i1 = i;
	} else {
	  *i0 = i;
	}
      }
    } 
  } else {    
    if (x >= xg[0]) {
      *i0 = 0;
      *i1 = 1;
    } else if (x <= xg[n-1]) {
      *i0 = n-1;
      *i1 = n-2;
      f = 0.0;
    } else {
      *i0 = 0;
      *i1 = n-1;
      while (*i1 - *i0 > 1) {
	i = (*i0 + *i1)/2;
	if (x >= xg[i]) {
	  *i1 = i;
	} else {
	  *i0 = i;
	}
      }
    } 
  }
  f = (x - xg[*i0])/(xg[*i1] - xg[*i0]);
  
  return f;
}

double interpol(double x, int n, double *xg, double *yg) {
  double f;
  int i0, i1;

  f = interpolx(x, n, xg, &i0, &i1);
  return (1.0-f)*yg[i0] + f*yg[i1];
}
 
double interpol2d(double x, double y, 
		  int nx, double *xg, 
		  int ny, double *yg,
		  double **z) {
  double fx, fy, z0, z1;
  int ix0, ix1, iy0, iy1;
  
  fx = interpolx(x, nx, xg, &ix0, &ix1);
  fy = interpolx(y, ny, yg, &iy0, &iy1);
  z0 = (1-fx)*z[ix0][iy0] + fx*z[ix1][iy0];
  z1 = (1-fx)*z[ix0][iy1] + fx*z[ix1][iy1];
  return (1-fy)*z0 + fy*z1;
}

double net_hc(HCF *hcf, double logxi, double logt) {
  double h, c;

  if (logxi < XMIN) logxi = XMIN;
  if (logxi > XMAX) logxi = XMAX;
  if (logt < TMIN) logt = TMIN;
  if (logt > TMAX) logt = TMAX;  

  h = interpol2d(logxi, logt, hcf->nx, hcf->xi, hcf->nt, hcf->te, hcf->htot);
  c = interpol2d(logxi, logt, hcf->nx, hcf->xi, hcf->nt, hcf->te, hcf->ctot);
  h = pow(10, h);
  c = pow(10, c);

  /*
  if (logt <= hcf->tmin) {
    double c0 = interpol2d(logxi, hcf->tmin, hcf->nx, hcf->xi, hcf->nt, hcf->te, hcf->ctot);
    double h0 = interpol2d(logxi, hcf->tmin, hcf->nx, hcf->xi, hcf->nt, hcf->te, hcf->htot);
    if (h0 < c0) {
      c0 = pow(10, c0);
      h0 = pow(10, h0);
      h += c0-h0;
    }
    if (h < c) return (h-c)*(logt-hcf->tmin);
  }
  */
  
  return h-c;
}

double mwgt(HCF *hcf, double logxi, double logt) {
  if (logxi < XMIN) logxi = XMIN;
  if (logxi > XMAX) logxi = XMAX;
  if (logt < TMIN) logt = TMIN;
  if (logt > TMAX) logt = TMAX;
  return interpol2d(logxi, logt, hcf->nx, hcf->xi, hcf->nt, hcf->te, hcf->mwgt);
}

double racc_line(HCF *hcf, double logxi, double logt) {
  if (logxi < XMIN) logxi = XMIN;
  if (logxi > XMAX) logxi = XMAX;
  if (logt < TMIN) logt = TMIN;
  if (logt > TMAX) logt = TMAX;
  double r=interpol2d(logxi, logt, hcf->nx, hcf->xi, hcf->nt, hcf->te, hcf->racc1);
  r = pow(10, r);
  return r;
}

double racc_abs(HCF *hcf, double logxi, double logt) {
  if (logxi < XMIN) logxi = XMIN;
  if (logxi > XMAX) logxi = XMAX;
  if (logt < TMIN) logt = TMIN;
  if (logt > TMAX) logt = TMAX;
  double r=interpol2d(logxi, logt, hcf->nx, hcf->xi, hcf->nt, hcf->te, hcf->racc2);
  r = pow(10, r);
  return r;
}

double racc_sct(HCF *hcf, double logxi, double logt) {
  if (logxi < XMIN) logxi = XMIN;
  if (logxi > XMAX) logxi = XMAX;
  if (logt < TMIN) logt = TMIN;
  if (logt > TMAX) logt = TMAX;
  double r=interpol2d(logxi, logt, hcf->nx, hcf->xi, hcf->nt, hcf->te, hcf->racc3);
  r = pow(10, r);
  return r;
}

double racc_tot(HCF *hcf, double logxi, double logt) {
  if (logxi < XMIN) logxi = XMIN;
  if (logxi > XMAX) logxi = XMAX;
  if (logt < TMIN) logt = TMIN;
  if (logt > TMAX) logt = TMAX;
  double r=interpol2d(logxi, logt, hcf->nx, hcf->xi, hcf->nt, hcf->te, hcf->racc);
  r = pow(10, r);
  return r;
}

double ion_frac(HCF *hcf, double logxi, double logt) {
  if (logxi < XMIN) logxi = XMIN;
  if (logxi > XMAX) logxi = XMAX;
  if (logt < TMIN) logt = TMIN;
  if (logt > TMAX) logt = TMAX;
  double r;
  r = interpol2d(logxi, logt, hcf->nx, hcf->xi, hcf->nt, hcf->te, hcf->fion);
  return r;
}
  
void set_esrc0(VH1State *vh1s) {
  double hdt = 0.5*vh1s->dt;
  for (int i = 0; i < vh1s->nrt; i++) {
    vh1s->esrc[i] = hdt*((vh1s->f0[i]+vh1s->g0[i]+vh1s->rf[i])*vh1s->u0[i] + 
			 (vh1s->f1[i]+vh1s->g1[i]+vh1s->rf[i])*vh1s->u[i]);
  }
}

void set_lumr(VH1State *vh1s) {
  int j, nr;
  double u, d, n, r, ru, t0, logt0, logx, fes;
  double nru2, n2u, r2p4, rd, hr, dlum, lum;
  HCF hc;

  if (vh1s->lum < 0) vh1s->lum = 0.0;
  for (int i = 0; i <= NGHOST ; i++) {
    vh1s->lumr[i] = vh1s->lum;
    vh1s->lumrs[i] = vh1s->lum;
    vh1s->xi[i] = 0.0;
    vh1s->hr[i] = 0.0;
  }
  u = vh1s->ulum/vh1s->uvol;
  int im = vh1s->nrt-1;
  logx = XMIN;
  logt0 = TMIN;
  for (int i=NGHOST; i < vh1s->nrt-1; i++) {
    d = vh1s->d[i];
    n = (d*vh1s->uden*vh1s->fhm)/MP;
    r = vh1s->rc[i];
    ru = r*vh1s->ulength; 
    nru2 = vh1s->ulum/(n*ru*ru);
    n2u = n*n/u;
    r2p4 = 4*PI*r*r;
    vh1s->xi[i] = vh1s->lumrs[i]*nru2;
    logx = vh1s->xi[i]>0?log10(vh1s->xi[i]):XMIN;
    t0 = qtemp(vh1s->hcf, vh1s->p[i]/vh1s->d[i], vh1s->uvel, logx);
    vh1s->te[i] = t0;
    logt0 = t0>0?log10(t0):TMIN;
    vh1s->fion[i] = ion_frac(vh1s->hcf, logx, logt0);
    hr = racc_abs(vh1s->hcf, logx, logt0)*n;
    hr *= (d*vh1s->uden)*C/u;
    dlum = r2p4*hr;
    rd = vh1s->lumrs[i]*0.05/dlum;
    nr = 1+(int)(vh1s->rd[i]/rd);
    if (nr < 25) nr = 25;
    if (nr > 250) nr = 250;
    rd = vh1s->rd[i]/nr;
    vh1s->hr[i] = racc_tot(vh1s->hcf, logx, logt0)*n;
    vh1s->hr[i] *= (d*vh1s->uden)*C/u;
    lum = vh1s->lumrs[i];
    vh1s->lumr[i] = lum;
    if (vh1s->lumrs[i] <= 0) {
      im = i;
      break;
    }
    dlum = lum*(1-exp(-dlum*rd/lum));
    for (j = 1; j < nr; j++) {
      lum -= dlum;
      if (lum <= 0) {
	break;
      }
      vh1s->lumr[i] += lum;
      vh1s->xi[i] = lum*nru2;    
      logx = vh1s->xi[i]>0?log10(vh1s->xi[i]):XMIN;
      hr = racc_abs(vh1s->hcf, logx, logt0)*n;
      hr *= (d*vh1s->uden)*C/u;
      dlum = lum*(1-exp(-r2p4*hr*rd/lum)); 
      hr = racc_tot(vh1s->hcf, logx, logt0)*n;
      hr *= (d*vh1s->uden)*C/u;
      vh1s->hr[i] += hr;
    }
    if (lum > dlum) {
      vh1s->lumrs[i+1] = lum-dlum;
    } else {
      vh1s->lumrs[i+1] = 0.0;
    }
    vh1s->lumr[i] /= nr;
    vh1s->hr[i] /= nr;
  }

  int i = im;  
  for (j = i+1; j < vh1s->nrt; j++) {
    vh1s->lumrs[j] = 0.0;
    vh1s->lumr[j] = 0.0;
    vh1s->hr[j] = 0.0;
    vh1s->xi[j] = vh1s->xi[i];
    t0 = qtemp(vh1s->hcf, vh1s->p[i]/vh1s->d[i], vh1s->uvel, logx);
    vh1s->te[i] = t0;
    logt0 = t0>0?log10(t0):TMIN;
    vh1s->fion[i] = ion_frac(vh1s->hcf, logx, logt0);
    /*
    vh1s->lumrs[j] = 0;
    vh1s->lumr[j] = 0;
    vh1s->xi[j] = 0;
    vh1s->hr[j] = 0;
    vh1s->fion[j] = 0;
    */
  }
  if (i <= NGHOST) {
    for (j = 0; j < i; j++) {
      vh1s->lumrs[j] = 0;
      vh1s->lumr[j] = 0;
      vh1s->xi[j] = 0;
      vh1s->hr[j] = 0;
      vh1s->fion[j] = 0;
    }
  } else {
    for (i = 0; i < NGHOST; i++) {
      vh1s->lumrs[i] = vh1s->lumrs[NGHOST];
      vh1s->lumr[i] = vh1s->lumr[NGHOST];
      vh1s->xi[i] = vh1s->xi[NGHOST];
      vh1s->hr[i] = vh1s->hr[NGHOST];
      vh1s->fion[i] = vh1s->fion[NGHOST];
    }
    j = vh1s->nrt-NGHOST-1;
    for (i = j+1; i < vh1s->nrt; i++) {
      vh1s->lumrs[i] = vh1s->lumrs[j];
      vh1s->lumr[i] = vh1s->lumr[j];
      vh1s->xi[i] = vh1s->xi[j];
      vh1s->hr[i] = vh1s->hr[j];
      vh1s->fion[i] = vh1s->fion[j];
    }
  }
}

void set_esrc1ev(VH1State *vh1s) {
  double u = vh1s->ulum / vh1s->uvol;
  double ue = vh1s->uenergy/vh1s->umass;
  double d, n, r, x, q0, q1, t, h, logt, y0, t0, logt0, qe, logx, ud;
  for (int i = 0; i < vh1s->nrt; i++) {
    d = vh1s->d[i];
    n = (d*vh1s->fhm*vh1s->uden)/MP;
    r = vh1s->rc[i];
    r *= vh1s->ulength;
    x = (vh1s->lumr[i]*vh1s->ulum)/(n*r*r);
    vh1s->xi[i] = x;
    logx = (x>0)?log10(x):XMIN;
    if (logx > XMAX) logx = XMAX;
    ud = n*n/(u*d);
    q0 = vh1s->p[i]/(vh1s->gamma1*d);
    q1 = q0;
    t = 0;
    h = vh1s->dt*0.25;
    while (t < vh1s->dt) {
      logt0 = interpol(logx, vh1s->hcf->nx, vh1s->hcf->xi, vh1s->hcf->ete);
      if (logt0 < vh1s->hcf->tmin) {
	logt0 = vh1s->hcf->tmin;
      }
      t0 = pow(10, logt0);
      y0 = mwgt(vh1s->hcf, logx, logt0);
      qe = ((KB*t0)/(y0*MP*vh1s->gamma1))/ue;
      if (q0 < qe) {
	q1 = qe;
	t = vh1s->dt;
	break;
      }
      int nit = 0;
      while (1) {
	nit++;
	double te = qtemp(vh1s->hcf, q1*vh1s->gamma1, vh1s->uvel, logx);
	double hc = ud*net_hc(vh1s->hcf, logx, log10(te));
	double dq = hc*(vh1s->dt-t);
	double q1s = dq*0.1;
	double dqe = qe-q1;
	if (((dq > 0) && dqe <= q1s) ||
	    ((dq < 0) && dqe >= q1s)) {
	  q1 = qe;
	  t = vh1s->dt;
	  break;
	} else {
	  q1s = 0.01*fabs(dqe);
	  if (fabs(dq) <= q1s) {
	    q1 += dq;
	    t = vh1s->dt;
	    break;
	  } else {
	    h = q1s/fabs(hc);
	    double tn = t+h;
	    if (tn >= vh1s->dt) {
	      q1 += hc*(vh1s->dt-t);
	      t = vh1s->dt;
	      break;
	    } else {
	      if (hc > 0) q1 += q1s;
	      else q1 -= q1s;
	      t = tn;
	    }
	  }
	}
	
	if (nit > 100000) {
	  printf("max iteration in hc inner loop: %lld %d %g %g %g %g %g %g %g\n",
		 vh1s->ncycle, i, q0, q1, te, t, vh1s->dt, q1s, hc);
	  exit(1);
	}
      }
    }  
    vh1s->esrc[i] = q1-q0;
    vh1s->te[i] = qtemp(vh1s->hcf, q1*vh1s->gamma1, vh1s->uvel, logx);
    logt = (vh1s->te[i] > 0)?log10(vh1s->te[i]):TMIN;
    vh1s->fion[i] = ion_frac(vh1s->hcf, logx, logt);
    vh1s->lumc[i] = 4*PI*vh1s->rc[i]*vh1s->rc[i]*vh1s->rd[i]*vh1s->d[i]*(q1-q0)/vh1s->dt;
  }
}

void set_esrc1eq(VH1State *vh1s) {
  int j, nt, m;
  double n, r, x, t, u, d, s;
  HCF hc;
  double hdt, tt, dtt, dt, logt0, t0, t1, q0, q1;
  double ud, ue, y0, y1, qt, dhc, dtp, logx, qp, tp, logtp;
  double dtt0;  
  if (vh1s->teq) {
    dtt0 = 0.5;
  } else {
    dtt0 = 0.05;
  }  
  u = vh1s->ulum/vh1s->uvol;
  ue = vh1s->uenergy/vh1s->umass;
  dt = vh1s->dt;
  for (int i = NGHOST-3; i < vh1s->nrt-3; i++) {
    d = vh1s->d[i];
    n = (d*vh1s->fhm*vh1s->uden)/MP;
    r = vh1s->rc[i];
    r *= vh1s->ulength; 
    x = (vh1s->lumr[i]*vh1s->ulum)/(n*r*r);
    logx = x>0?log10(x):XMIN;
    if (logx < XMIN) logx = XMIN;
    if (logx > XMAX) logx = XMAX;
    q0 = q1 = vh1s->p[i]/(vh1s->gamma1*vh1s->d[i]);
    if (vh1s->teq == 2) {
      logt0 = interpol(logx, vh1s->hcf->nx, vh1s->hcf->xi, vh1s->hcf->ete);
      if (logt0 < vh1s->hcf->tmin) {
	logt0 = vh1s->hcf->tmin;
      }
      t0 = pow(10, logt0);
      y0 = mwgt(vh1s->hcf, logx, logt0);
      q0 = ((KB*t0)/(y0*MP*vh1s->gamma1))/ue;
      dhc = q0-q1;
      vh1s->esrc[i] = dhc;
      t0 = qtemp(vh1s->hcf, q0*vh1s->gamma1, vh1s->uvel, logx);
      logt0 = log10(t0);
      vh1s->te[i] = t0;
      vh1s->fion[i] = ion_frac(vh1s->hcf, logx, logt0);
      vh1s->lumc[i] = 4*PI*vh1s->rc[i]*vh1s->rc[i]*vh1s->rd[i]*vh1s->d[i]*(q0-q1)/vh1s->dt;
      continue;
    }
    logt0 = vh1s->te[i] > 0? log10(vh1s->te[i]) : TMIN;
    t0 = qtemp(vh1s->hcf, q0*vh1s->gamma1, vh1s->uvel, logx);
    
    if (1+x==x) {      
      vh1s->esrc[i] = 0.0;
      vh1s->te[i] = t0;
      continue;
    }

    if (logx <= XMIN && vh1s->te[i] <= TMIN) {
      vh1s->esrc[i] = 0.0;
      continue;
    }

    t1 = t0;
    logt0 = log10(t0);
    vh1s->xi[i] = x;
    ud = n*n/(u*d);
    tt = 0.0;
    y0 = net_hc(vh1s->hcf, logx, logt0);
    double ys = y0;
    double y1s;
    if (y0 != 0) {
      m = 0;
      if (y0 > 0) {
	dtt = dtt0;
      } else {
	dtt = -dtt0;
      }
      y0 = 1.0/(y0*ud);
      dtp = 0;
      while (1) {
	if ((logt0 < TMIN && dtt < 0) || (logt0 > TMAX) && ( dtt > 0)) {
	  printf("too small or too large T: %lld %d %g %g %g %g %g %g %g %g %g %g\n", 
		 vh1s->ncycle, i, q1, q0, t1, t0, dtt, x, logx, vh1s->lumr[i], vh1s->lumrs[i], n);
	  break;
	}
	if (m > vh1s->miter) {
	  printf("too many iterations in esrc1: %lld %d %g %g %g %g %g %g %g %g\n",
		 vh1s->ncycle, i, q1, q0, t1, t0, logx, tt, dt, dtt);
	  break;
	}
	qp = q0*exp(dtt);
	tp = qtemp(vh1s->hcf, qp*vh1s->gamma1, vh1s->uvel, logx);
	logtp = log10(tp);
	y1 = net_hc(vh1s->hcf, logx, logtp);
	y1s = y1;
	if (y1 != 0) {
	  y1 = 1.0/(y1*ud);
	  if (y1*dtt < 0) {
	    dtt /= 2.0;
	    m++;
	    continue;
	  }
	  if (vh1s->teq) {
	    if (fabs(dtt) < 1e-3) break;
	  } else {
	    double ay = 0.5*(y0+y1);
	    double dy = fabs(y1-y0)/ay;
	    if (dy < 0.05) {
	      dtp = ay*(qp-q0);
	      if (fabs(dtt) < 1e-3 && dtp/dt < 1e-3) break;
	      if (tt+dtp > dt) {	      
		dtt /= 2.0;
		m++;
		continue;
	      } else if (dy < 0.01 && fabs(dtt)<1e-3 && (dt-tt)>5*dtp) {
		dtt *= 2;
	      }
	    } else {
	      dtt /= 2.0;
	      m++;
	      continue;
	    }
	  }
	} else {
	  if (vh1s->teq) {
	    break;
	  } else {
	    dtt /= 2.0;
	    m++;
	    continue;
	  }
	}
	y0 = y1;
	q0 = qp;
	t0 = tp;
	logt0 = logtp;
	tt += dtp;
      }
    }
    dhc = q0-q1;
    vh1s->esrc[i] = dhc;
    vh1s->te[i] = t0;
    logt0 = (t0 > 0)?log10(t0):TMIN;
    vh1s->fion[i] = ion_frac(vh1s->hcf, logx, logt0);
    vh1s->lumc[i] = 4*PI*vh1s->rc[i]*vh1s->rc[i]*vh1s->rd[i]*vh1s->d[i]*(q0-q1)/vh1s->dt;
  }  
  /*
  for (i = NGHOST+1; i < vh1s->nrt; i++) {
    dhc = -vh1s->esrc[i-1];
    if (dhc > 0) { 
      dhc *= 4*PI*vh1s->rc[i]*vh1s->rc[i]*vh1s->rd[i]*vh1s->d[i];    
      dhc /= vh1s->dt;
    } else {
      dhc = 0;
    }
    vh1s->lumc[i] = vh1s->lumc[i-1] + dhc;
  }
  */
}

void mdot_stat(VH1State *vh1s) {
  int i, j = NGHOST;
  vh1s->mdot =  -4*PI*vh1s->rc[j]*vh1s->rc[j]*vh1s->d[j]*vh1s->u[j];
  j = vh1s->nr + NGHOST-1;
  vh1s->mdoto =  -4*PI*vh1s->rc[j]*vh1s->rc[j]*vh1s->d[j]*vh1s->u[j];
  //if (vh1s->mdot > 0) {
  //  vh1s->mb += vh1s->mdot * vh1s->dt;
  //}
  double c2 = pow(C/vh1s->uvel, 2);

  if (vh1s->mdot > 0) {
    if (vh1s->reffs >= 0.0) {
      double a = pow(0.9663-0.9292*vh1s->reffs, -0.5639);
      double b = pow(4.6270-4.4450*vh1s->reffs, -0.5524);
      double c = pow(827.3-718.1*vh1s->reffs, -0.706);
      double rx = vh1s->mdot/(16.0*vh1s->medd);      
      double rc = vh1s->reffc;
      if (rc <= 0) rc = 1.0;
      if (vh1s->refft > 0 && rc > vh1s->reffd) {
	rc = vh1s->reffd + (rc-vh1s->reffd)*(1-exp(-vh1s->time/vh1s->refft));
      }
      if (rx < 1e-10) rx = 1e-10;
      rx = 1/rx;
      double lam = (a*(0.985/(rx+b) + 0.015/(rx+c)))*rc;
      vh1s->reff = (rx/16.0)*lam;
    } else if (vh1s->reffa > 0 && vh1s->reffb > 0) {
      double md = vh1s->mdot/vh1s->medd;
      double rc = vh1s->reffc;
      if (rc <= 0) rc = 1.0;
      if (vh1s->refft > 0 && rc > vh1s->reffd) {
	rc = vh1s->reffd + (rc-vh1s->reffd)*(1-exp(-vh1s->time/vh1s->refft));
      }
      md /= rc;
      vh1s->reff = vh1s->reffa/(vh1s->reffb + vh1s->reffa*md);
    }
  }
  if (vh1s->mdot > 0) {
    vh1s->lumi = vh1s->reff * vh1s->mdot * c2;
  } else {
    vh1s->lumi = 0.0;
  }
  if (vh1s->lumts <= 0) {
    vh1s->lum = vh1s->lumi;
  } else {
    double t1 = vh1s->time/vh1s->lumts;
    double t21 = vh1s->dt/vh1s->lumts;
    double t2 = t1+t21;
    if (t21 <= 0) {
      vh1s->lum = vh1s->lumi;
    } else {
      t21 = exp(-t21);
      t1 = exp(-t1);
      t2 = exp(-t2);
      vh1s->lum = (t21*(1-t1)*vh1s->lum + (1-t21)*vh1s->lumi)/(1-t2);
    }
  }
  vh1s->emb = vh1s->mb * (1 - vh1s->lum/vh1s->lumedd);
}

double nc_scale(double m) {
  return 2e10/m;
}

double mdot_scale(double n, double m) {
  double x = n/nc_scale(m);
  return pow(x/(1+x),0.5);
}

double ri_scale(double n, double m) {
  double md = mdot_scale(n, m);
  return 338*pow(md*0.01*m/(n*n), 1.0/3.0);
}

double rb_scale(double m) {
  return 6.5e-3*0.01*m;
}

double tc_scale(double n, double m) {
  return ri_scale(n,m)/3.3e-5;
}

double wall_time() {
  clock_t t;
  t = clock();
  double r = ((double)t)/CLOCKS_PER_SEC;
  return r;
}
  
void evolve_vh1state(VH1State *vh1s) {
  FILE *fd, *f;
  char fname[1000];
  int i, j;
  double dm, db, pb, ub, db0, pb0, ub0;
  double lastmd, lastt, md1, md2, tcycle;
  long long lastcycle;

  fd = fopen(vh1s->fmdot, "w");
  if (fd == NULL) {
    printf("cannt open file %s\n", vh1s->fmdot);
    exit(1);
  }
  j = NGHOST;
  vh1s->dt = 0.0;
  vh1s->mtoti = 0.0;
  vh1s->mtoto = 0.0;
  mdot_stat(vh1s);
  db0 = vh1s->dib/vh1s->uden;
  pb0 = vh1s->pi/(vh1s->uvel*vh1s->uvel*vh1s->uden);
  ub0 = vh1s->vi/vh1s->uvel;
  db = db0;
  pb = pb0;
  ub = ub0;
  lastcycle = 0;
  lastmd = 0.0;
  lastt = vh1s->time;
  double minmd = vh1s->mdot;
  double maxmd = vh1s->mdot;
  double lastwt = wall_time();
  double wt0 = lastwt;
  tcycle = vh1s->tcycle;
  tcycle *= YEAR/vh1s->utime;
  md1 = 0.0;
  md2 = 0.0;  
  int stopsol = 0;
  long long ncycle0 = vh1s->ncycle;
  double time0 = vh1s->time;
  double tw1, tw2, tw3, tw4, tw5, tw6, tw7, tw8;
  tw1 = 0.0;
  tw2 = 0.0;
  tw3 = 0.0;
  tw4 = 0.0;
  tw5 = 0.0;
  tw6 = 0.0;
  double nlumupd = 0;
  double lumupd = 1E31;
  long long updcycle = ncycle0;
  int ilumupd = 1;
  long long endcycle = 0;
  if (vh1s->endcycle > 0) {
    endcycle += ncycle0 + vh1s->endcycle;
  }
  while (vh1s->time < vh1s->endtime &&
	 (endcycle <= 0 || vh1s->ncycle < endcycle)) {
    if (vh1s->ncycle%vh1s->ncheck == 0) {
      FILE *fs = fopen("stop_solution", "r");
      if (fs != NULL) {
	printf("stop printing solutions\n");
	stopsol = 1;
	system("rm stop_solution");
	fclose(fs);
      }
      fs = fopen("start_solution", "r");
      if (fs != NULL) {
	printf("start printing solutions\n");
	stopsol = 0;
	system("rm start_solution");
	fclose(fs);
      }
    }
    if (vh1s->mgasm > 0) {
      vh1s->mgas = 0.0;
      for (i = NGHOST; i < vh1s->nr+NGHOST; i++) {
	vh1s->mgas += 4*PI*vh1s->rc[i]*vh1s->rc[i]*vh1s->rd[i]*vh1s->d[i];
      }
      dm = vh1s->mgas - vh1s->mgas0;
      dm = exp(-dm/vh1s->mgasm);
      db = db0 * dm;
      pb = pb0 * dm;
      ub = ub0;
      i = 1;
      setbd_(&i, &db, &pb, &ub);
    }
    double w0 = wall_time();
    if (ilumupd) {
      set_lumr(vh1s);
    }
    double w1 = wall_time();
    if (ilumupd) {
      set_gsrc0(vh1s);
    }
    set_dt(vh1s);
    setgsrc_(&vh1s->nrt, vh1s->gsrc);
    double w2 = wall_time();
    sweepns_();
    oldsta_(&vh1s->nrt, vh1s->d0, vh1s->p0, vh1s->u0);    
    newsta_(&vh1s->nrt, vh1s->d, vh1s->p, vh1s->u, vh1s->nrc);
    double w3 = wall_time();
    double x = log10(Max(XMIN,vh1s->xi[j]));
    double te0 = qtemp(vh1s->hcf, vh1s->p[j]/vh1s->d[j], vh1s->uvel, x);
    set_gsrc1(vh1s);
    setgsrc_(&vh1s->nrt, vh1s->gsrc);
    updu_();
    //oldsta_(&vh1s->nrt, vh1s->d0, vh1s->p0, vh1s->u0);
    //newsta_(&vh1s->nrt, vh1s->d, vh1s->p, vh1s->u, vh1s->nrc);
    //printf("s1: %d %g %g %g %g\n", vh1s->ncycle, vh1s->p[j], vh1s->d[j], vh1s->te[j], te0);
    //set_esrc0(vh1s);
    //setesrc_(&vh1s->nrt, vh1s->esrc);    
    //upde_();
    //newsta_(&vh1s->nrt, vh1s->d, vh1s->p, vh1s->u, vh1s->nrc);
    //printf("s1: %d %g %g %g %g\n", vh1s->ncycle, vh1s->p[j], vh1s->d[j], vh1s->te[j], te0);
    //double w5 = wall_time();
    remap_();
    newsta_(&vh1s->nrt, vh1s->d, vh1s->p, vh1s->u, vh1s->nrc);
    double w4 = wall_time();
    double te1 = qtemp(vh1s->hcf, vh1s->p[j]/vh1s->d[j], vh1s->uvel, x);
    double te2 = te1;
    if (ilumupd) {
      if (vh1s->teq <= 2) {
	set_esrc1eq(vh1s);
      } else {
	set_esrc1ev(vh1s);
      }
      setesrc_(&vh1s->nrt, vh1s->esrc);
      updq_();
      newsta_(&vh1s->nrt, vh1s->d, vh1s->p, vh1s->u, vh1s->nrc);
      te2 = qtemp(vh1s->hcf, vh1s->p[j]/vh1s->d[j], vh1s->uvel, x);
    }
    double w5 = wall_time(); 
    if (vh1s->mdot > 0) {
      vh1s->mtoti += vh1s->dt*vh1s->mdot;
    } else {
      vh1s->mtoto -= vh1s->dt*vh1s->mdot;
    }
    if (vh1s->mdot > maxmd) maxmd = vh1s->mdot;
    if (vh1s->mdot < minmd) minmd = vh1s->mdot;
    mdot_stat(vh1s);
    double mdr = 0.5*(vh1s->lum+lumupd);
    double mdn = vh1s->lum;
    if (vh1s->ncycle-updcycle >= vh1s->lum_updc ||
	fabs(mdn-lumupd) > vh1s->lum_updm*mdr) {
      ilumupd = 1;
      lumupd = mdn;
      updcycle = vh1s->ncycle;
      nlumupd += 1;
    } else {
      ilumupd = 0;
    }
    double wt = wall_time();
    tw1 += w1-w0;
    tw2 += w2-w1;
    tw3 += w3-w2;
    tw4 += w4-w3;
    tw5 += w5-w4;
    tw6 += wt-w5;
    double dt0 = vh1s->dt;
    long long thiscycle = vh1s->ncycle;
    double thistime = vh1s->time;
    int smalldt = 0;
    if (smalldt || (thiscycle == ncycle0+1) ||
	(thiscycle-lastcycle == vh1s->nmdotmax) ||
	(thistime-lastt >= vh1s->tmdotmax*tcycle) ||
	((thiscycle-lastcycle >= vh1s->nmdotmin) &&
	 (thistime-lastt >= vh1s->tmdotmin*tcycle) &&
	 fabs(lastmd-vh1s->mdot) > vh1s->vmdot*fabs(vh1s->mdot)) ||
	(vh1s->wmdot > 0 && wt-lastwt >= vh1s->wmdot)) {
      fprintf(fd, "%10d %15.8E %15.8E %5d %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %9.3E %9.3E %9.3E %9.3E %9.3E %9.3E %9.3E %9.3E %9.3E %d\n", 
	      thiscycle, thistime, dt0, j, vh1s->rc[j], 
	      vh1s->d[j], vh1s->p[j], vh1s->u[j], 
	      vh1s->mdot, vh1s->mtoti/thistime,
	      vh1s->mtoto/thistime, maxmd, minmd,
	      vh1s->mb, vh1s->emb, vh1s->lumedd,
	      vh1s->lumi, vh1s->lum, vh1s->reff,
	      wt-wt0, tw1, tw2, tw3, tw4, tw5, tw6, wt-lastwt,
	      nlumupd/(thiscycle-lastcycle), (int)(0.1+nlumupd));
      fflush(fd);
      if (smalldt || !stopsol) {
	sprintf(fname, "%s%d.txt", vh1s->fsolution, thiscycle);
	f = fopen(fname, "w");
	if (f == NULL) {
	  printf("cannot open file %s\n", fname);
	  exit(1);
	}	
	fprintf(f, "# %10d %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E\n",
		thiscycle, thistime, dt0,
		vh1s->mdot, vh1s->mtoti, vh1s->mtoto,
		vh1s->mb, vh1s->emb, vh1s->lumedd,
		vh1s->lumi, vh1s->lum, vh1s->reff);
	for (i = 0; i < vh1s->nrt; i++) {
	  double logxi = vh1s->xi[i]>0?log10(vh1s->xi[i]):XMIN;
	  double logti = vh1s->te[i]>0?log10(vh1s->te[i]):TMIN;
	  double mmw = mwgt(vh1s->hcf, logxi, logti);
	  fprintf(f, "%10d %15.8E %5d %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E\n", 
		  thiscycle, thistime, i, vh1s->rc[i], 
		  vh1s->p[i], vh1s->d[i], vh1s->u[i], vh1s->te[i],
		  vh1s->f0[i], vh1s->g0[i], vh1s->rf[i], 
		  vh1s->lumr[i], vh1s->xi[i], vh1s->fion[i],
		  vh1s->hr[i], mmw, vh1s->lumc[i]);
	}
	fclose(f);
      }
      lastcycle = thiscycle;
      lastmd = vh1s->mdot;
      lastt = thistime;
      lastwt = wt;
      tw1 = 0.0;
      tw2 = 0.0;
      tw3 = 0.0;
      tw4 = 0.0;
      tw5 = 0.0;
      tw6 = 0.0;
      nlumupd = 0;
    }
    md2 = md1;
    md1 = vh1s->mdot;
    if (smalldt) {
      printf("timestep has become too small: %d %g %g %g %g\n", thiscycle, thistime, dt0, vh1s->dt);
      break;
    }
  }
}

