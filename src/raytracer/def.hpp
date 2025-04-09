#ifndef _DEF_H
#define _DEF_H

#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
#include <cmath>
#include <iostream>
using namespace std;

#define imax 400

const long double Pi = acos(-1.0);

long double xobs, yobs;
long double epsi3, a13, a22, a52;
long double spin;
long double iobs_deg;
int phicount;

struct SurfacePoint {
  long double x;
  long double y;
  long double u0, u1, u2, u3;            
};

struct RayHit {
  long double cosem;
  long double gfactor;
  long double r;
};

/*-----------------------------------------------------------*/

void raytrace(long double xobs, long double yobs, long double iobs,
              long double rin, long double disk_length_combined,
              RayHit &hit, int &stop_integration,
              const SurfacePoint *diskdata, const size_t ddsize);
void diffeqs(long double b, long double vars[], long double diffs[]);
void redshift(long double r, long double ktkp, long double &gg);
// void redshift_polish_doughnut(long double r, long double th, long double l
// ,long double ktkp, long double& gg);
void intersection(long double x_1, long double y_1, long double z_1,
                  long double x_2, long double y_2, long double z_2,
                  long double x_d[]);
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
