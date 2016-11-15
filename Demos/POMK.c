#include "auto_f2c.h"
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/*  mk model :       Modified Klausmeier                                  */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int func (integer ndim, const doublereal *u, const integer *icp,
          const doublereal *par, integer ijac,
          doublereal *f, doublereal *dfdu, doublereal *dfdp)
{
  /* Local variables */
  doublereal B,W,Bx,Wx, yy,zz,L;
  doublereal a,d,m;

  a    = par[0];
  d    = par[1];
  m    = par[2];

  yy     = par[18];
  zz 	 = par[19];
  L  	 = par[20];

  B      = u[0];
  W      = u[1];
  Bx     = u[2];  
  Wx     = u[3];

  /* Define the first derivatives */
  f[0] = L * Bx;
  f[1] = L * Wx;
    
  /* Define the two equations (divided by the respective diffusion parameter) */
  f[2] = -L * ( W*B*B - m*B + zz*Bx) / 1;
  f[3] = -L * ( a - W - W*B*B  + (yy+zz)*Wx) / d;

  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int stpnt (integer ndim, doublereal t,
           doublereal *u, doublereal *par)
{
  doublereal B0,W0;
  doublereal a,d,m,L;

  /* Initialize the equation parameters */
  a = 1.1;
  d = 8;
  m = 0.5;
  L = 1.0;
    
  par[0] = a;
  par[1] = d;
  par[2] = m;
 
  /* Initilize dummy variables */
  par[18] = (doublereal).00000000000000000009;	/* yy  */
  par[19] = (doublereal).00000000000000000010;	/* zz  */
  par[20] = L;
    
  /* Solution of quadratic equation */
    W0  = a/2 - 0.5 * sqrt(a*a-4*m*m);
    B0  = m / W0;
 
  /* Initialize the solution */
  u[0] = B0;
  u[1] = W0;
  u[2] = 0.0;
  u[3] = 0.0;
 
  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int bcnd (integer ndim, const doublereal *par, const integer *icp,
          integer nbc, const doublereal *u0, const doublereal *u1, integer ijac,
          doublereal *fb, doublereal *dbc)
{
    fb[0] = u1[2];   // Bx_right   = 0
    fb[1] = u0[2];   // Bx_left    = 0
    fb[2] = u1[3];   // Wx_right   = 0
    fb[3] = u0[3];   // Wx_left    = 0

  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int icnd (integer ndim, const doublereal *par, const integer *icp,
          integer nint, const doublereal *u, const doublereal *uold,
          const doublereal *udot, const doublereal *upold, integer ijac,
          doublereal *fi, doublereal *dint)
{
  fi[0] = upold[0]*(u[0]-uold[0]);
  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int fopt (integer ndim, const doublereal *u, const integer *icp,
          const doublereal *par, integer ijac,
          doublereal *fs, doublereal *dfdu, doublereal *dfdp)
{
  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int pvls (integer ndim, const doublereal *u,
          doublereal *par)
{
  integer tmp;
  extern doublereal getp();
  extern doublereal getpuwe();

  if(par[10]>0)
  	par[11]=6.2832/par[10];
  else
	par[11]=0;
  par[12]=u[0]*u[0];

  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
