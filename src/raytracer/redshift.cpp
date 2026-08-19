#include <cmath>

#include "redshift.hpp"
#include "metric.hpp"

void redshift(Real r, Real ktkp, Real& gg) {
  /* I have to write this function new for the case of Polish doughnut model.*/
  /* Using equation 6.22 Cosimo's Book, and Ut and OMEGA from E.5 and E.6
   * respectively */
  Real Omega;
  Real uet;
  Real met[4][4], met_rder[4][4];
  Real th = Pi / 2;

  metric(r, th, met);

  metric_rderivatives(r, th, met_rder);

  Omega = (-met_rder[0][3] + std::sqrt((met_rder[0][3] * met_rder[0][3]) -
                                  (met_rder[0][0] * met_rder[3][3]))) /
          (met_rder[3][3]);

  // Omega  = (-rderg03 + sqrt(rderg03*rderg03 - rderg00*rderg33))/rderg33;

  // uet = sqrt(-g00 - 2.*g03*Omega - g33*Omega*Omega);

  uet =
      std::sqrt(-met[0][0] - (2 * met[0][3] * Omega) - (met[3][3] * Omega * Omega));

  gg = uet / (1. - ktkp * Omega);
}
