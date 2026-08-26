#include <cmath>

#include "redshift.hpp"
#include "metric.hpp"

void redshift_plunge(const Real& isco, const Real& r, Real* kvec, Real& gg){
  //Cardenas-Avendano2020
  Real ut_num = CUBE(spin)*r+SQR(spin)*std::sqrt(isco)*(-2*r+2*isco+r*isco)+spin*(CUBE(r)-2*SQR(isco))+CUBE(r)*std::sqrt(isco)*(isco-2.0);
  Real ut_den = r*(SQR(spin)+r*(r-2.0))*std::sqrt(2*spin*std::sqrt(CUBE(isco))-3*SQR(isco)+CUBE(isco));
  Real ut = ut_num/ut_den;

  Real ur_num_rad_num = r*SQR(r-isco)*(SQR(spin)*(-(r+2*isco))+4*spin*std::sqrt(isco)*(r+isco)+isco*(r*(isco-4)-2*isco));
  Real ur_num_rad_den = std::sqrt(CUBE(isco))*(SQR(spin)+(r-2)*r)*(2*spin+(isco-3)*std::sqrt(isco));
  Real ur_num_rad = - ur_num_rad_num/ur_num_rad_den;
  Real ur = -std::sqrt(SQR(spin)+(r-2)*r)/SQR(r)*std::sqrt(ur_num_rad);

  Real up_num = SQR(spin)*r+2*spin*(std::sqrt(CUBE(isco))-r*std::sqrt(isco))+SQR(isco)*(r-2.0);
  Real up_den = r*(SQR(spin)+r*(r-2))*std::sqrt(2*spin*std::sqrt(CUBE(isco))-3*SQR(isco)+CUBE(isco));
  Real up = up_num/up_den;

  gg = 1.0/(ut+up*kvec[3]/kvec[0]+ur*kvec[1]/kvec[0]);
}

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
