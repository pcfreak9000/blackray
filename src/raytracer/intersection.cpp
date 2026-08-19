#include "intersection.hpp"

#include <cmath>

void intersection(Real x_1, Real y_1, Real z_1,
                  Real x_2, Real y_2, Real z_2,
                  Real x_d[]) {
  Real r1, h1, r2, h2;
  Real pp;
  Real r3, r4, c1;

  r1 = x_1 * std::sin(y_1);
  h1 = x_1 * std::cos(y_1);

  r2 = x_2 * std::sin(y_2);
  h2 = -x_2 * std::cos(y_2);

  pp = std::cos(z_2 - z_1);

  r3 = std::sqrt(r1 * r1 + r2 * r2 - 2 * r1 * r2 * pp);

  r4 = h1 * r3 / (h1 + h2);

  c1 = (r1 - r2 * pp) / r3;

  x_d[0] = 0.;
  x_d[1] = std::sqrt(r1 * r1 + r4 * r4);

  if (x_d[1] > r1 && x_d[1] > r2) {
    x_d[1] = 0.5 * (r1 + r2);
  } else if (x_d[1] < r1 && x_d[1] < r2) {
    x_d[1] = 0.5 * (r1 + r2);
  } else { // printf("x_d[1] = %Le, r1 = %Le, r2 = %Le\n",x_d[1],r1,r2);
    x_d[2] = 0.5 * Pi;
  }

  x_d[3] = std::asin(r4 * std::sqrt(1. - c1 * c1) / x_d[1]);

  if (z_2 > z_1) {
    x_d[3] = z_1 + x_d[3];
  } else {
    x_d[3] = z_1 - x_d[3];
  }
}
