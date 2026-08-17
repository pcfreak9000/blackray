#pragma once

//Where is this used?
//#include <algorithm>

#define DEBUG_DIV 1.0
#define RING_DIV 15
#define MAX_ITER 10000
//Warn about approach to max iterations:
//#define ITER_WARN
//Warn about strongly deviating 4-velocity norms from the disk at emission points:
//#define DEBUG_FVEL_NORM
//Warn about strongly deviating 4-momentum norms of photons at emission points:
//#define DEBUG_FMOM_NORM
//Warn about and fix cosem>1.0.
#define DEBUG_COSEM
#define RESTRICT_DEBUGFILE_CRIT 0
#define DEBUGFILE_OUT_DIV 1

#define imax 400

using Real = long double;
#define NO_INTERSECT -1
#define INTERSECT 0
#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

const Real Pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282L;


extern Real epsi3, a13, a22, a52;
extern Real spin;
extern Real iobs_deg;
extern int phicount;

struct SurfacePoint {
  Real x;
  Real y;
  Real u0, u1, u2, u3;
  Real density;
};

struct SurfaceElement {
  SurfacePoint *sp0;
  SurfacePoint *sp1;
  int index;
};

struct RayHit {
  Real cosem;
  Real gfactor;
  Real r;
  Real hc;
};
