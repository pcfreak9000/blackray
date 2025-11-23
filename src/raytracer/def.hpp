#ifndef _DEF_H
#define _DEF_H

#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
using namespace std;

#define imax 400

const long double Pi = acos(-1.0);


long double epsi3, a13, a22, a52;
long double spin;
long double iobs_deg;
int phicount;

struct SurfacePoint {
  long double x;
  long double y;
  long double u0, u1, u2, u3;
  long double density;
};

struct SurfaceElement {
  SurfacePoint *sp0;
  SurfacePoint *sp1;
  int index;
};

struct RayHit {
  long double cosem;
  long double gfactor;
  long double r;
  long double hc;
};

using Real = long double;
#define Q1 0
#define Q2 1
#define Q3 2
#define Q4 3
#define NO_INTERSECT -1
#define INTERSECT 0
#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))


#define DEBUG_DIV 1.0
#define RING_DIV 15
#define MAX_ITER 6000
//Warn about approach to max iterations:
//#define ITER_WARN
//Warn about strongly deviating 4-velocity norms from the disk at emission points:
//#define DEBUG_FVEL_NORM
//Warn about strongly deviating 4-momentum norms of photons at emission points:
//#define DEBUG_FMOM_NORM
//Warn about and fix cosem>1.0.
//#define DEBUG_COSEM

class QuadTree {
public:
  QuadTree(Real x, Real y, Real width, Real height);
  Real check_intersect(Real x1, Real y1, Real x2, Real y2,
      SurfaceElement ** intersect);
  void put_element(SurfaceElement *element);
  void validate();
  size_t size();
private:
  std::vector<QuadTree*> subtrees;
  std::vector<SurfaceElement*> myelements;
  Real x, y, width, height;
  size_t max_elements;
  bool is_leaf;
  size_t level;
//  bool is_root;

  void subdivide();
  bool fits(SurfaceElement *element);
  bool overlaps(Real x0, Real y0, Real x1, Real y1);
  bool fully_inside(Real x0, Real y0, Real x1, Real y1);
};
/*-----------------------------------------------------------*/
Real checkIntersect(long double x1, long double y1, long double x2,
    long double y2, long double x3, long double y3, long double x4,
    long double y4);
void raytrace(long double xobs, long double yobs, long double iobs,
    long double rin, long double disk_length_combined, RayHit &hit,
    int &stop_integration, SurfacePoint **diskdata, const size_t ddsize, QuadTree* tree, Real checkr);
void diffeqs(long double b, long double vars[], long double diffs[]);
void redshift(long double r, long double ktkp, long double &gg);
// void redshift_polish_doughnut(long double r, long double th, long double l
// ,long double ktkp, long double& gg);
void intersection(long double x_1, long double y_1, long double z_1,
    long double x_2, long double y_2, long double z_2, long double x_d[]);
void metric(long double z1, long double z2, long double mn[][4]);
void metric_rderivatives(long double z1, long double z2, long double dmn[][4]);
void find_isco(long double z1, long double &isco);
// void polish_doughnut(long double r, long double theta, long double phi ,long
// double angm_disk, long double spin, long double& w_current); void
// emission_angle_rth(long double r, long double th, long double l, long double
// a, long double kr, long double kphi ,long double kth , long double lambda ,
// long double& em_angle);

#include "diffeqs.cpp"
#include "find_isco.cpp"
#include "intersection.cpp"
#include "metric.cpp"
#include "metric_rderivatives.cpp"
#include "raytracingnew.cpp"
#include "redshift.cpp"

#endif
