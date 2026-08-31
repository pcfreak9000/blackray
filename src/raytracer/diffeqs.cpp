#include "diffeqs.hpp"
#include <cmath>

void diffeqs(const Real& b, const Real *const vars, Real diffs[]) {
  const Real& r = vars[0];
  const Real& th = vars[1];

  const Real t1 = std::cos(th);
  const Real t2 = r * r; //r^2
  //t3 = pow(t2, 0.2e1); //r^4
  const Real t3 = SQR(t2);
  Real t4 = t2 * t3; //r^6
  Real t5 = r * t3; //r^5
  Real t6 = r * t2; //r^3
  const Real t7 = spin * spin;
  //t8 = t7 * pow(t1, 0.2e1);
  const Real t8 = t7 * SQR(t1);
  const Real t9 = (t8 + t2) * r + epsi3;
  const Real t10 = a22 + t2;
  const Real t11 = std::sin(th);
  const Real t12 = t7 + t2;
  //t13 = pow(t10, 0.2e1);
  const Real t13 = SQR(t10);
  //t14 = pow(t11, 0.2e1);
  const Real t14 = SQR(t11);
  const Real t15 = t7 * t13 * t14;
  Real t17 = a13 + t6;
  const Real t18 = t7 * r;
  Real t19 = t18 * t10 * t14;
  const Real t20 = t12 * t17;
  Real t21 = t20 - t19;
  const Real t22 = 2.0 * r;
  const Real t23 = -t22 + t12;
  //t17 = pow(t12, 0.2e1) * pow(t17, 0.2e1) - t7 * t4 * t23 * t14;
  t17 = SQR(t12) * SQR(t17) - t7 * t4 * t23 * t14;
  const Real t24 = a22 + t7;
  const Real t25 = r * a22;
  const Real t26 = t7 * a13;
  t4 = 2.0 * t4 + t2 * (a13 * t24 + (a22 * t7 + (a13 + t25) * r) * r) +
       t26 * a22;
  const Real t27 = a52 + t2;
  const Real t28 = t23 * t3 - t15;
  const Real t29 = t8 / 0.3e1 + t2;
  const Real t30 = 0.2e1 / 0.3e1 * t7;
  Real t31 = 0.3e1 / 0.5e1 * t7;
  const Real t32 = (t31 + t2) * r + 0.2e1 / 0.5e1 * a13;
  t31 = r * t32 - t31 * (a22 / 0.3e1 + t2) * t14;
  const Real t33 = 0.1e1 / t21;
  //t35 = pow(t33, 0.2e1);
  const Real t35 = SQR(t33);
  const Real t36 = t33 * t35;
  const Real t37 = 0.3e1 * t29;
  const Real t38 = t35 * t14;
  const Real t39 = t7 * a52;
  const Real t40 = a52 + t7;
  const Real t42 = 0.1e1 / r;
  const Real t43 = t9 * t42;
  const Real t44 =
      -2.0 * t11 * t35 * t1 * (-t43 * t17 + (t5 * t9 * t23 + t17) * t14 * t7) +
      4.0 * t1 * t7 * t10 * t11 * t14 * t9 * t17 * t36;
  t21 = 0.1e1 / t21;
  const Real t45 = 0.1e1 / t9;
  const Real t46 = 0.1e1 / t27;
  const Real t47 = 0.1e1 / t23;
  //t48 = pow(t21, 0.2e1);
  const Real t48 = SQR(t21);
  t19 = 2.0 * t4 * t1 * spin * t11 *
        (t18 * t14 *
             (t19 + a22 * epsi3 -
              ((-(a22 - t7) * r + a13 - epsi3) * r - t8 * t10) * r - t26) +
         t20 * t9) *
        t21 * t48;
  t21 = 2.0 * t7 * t1;
  t5 = r * t9 * (-t12 * t3 + 2.0 * t5 + t15) * t48;
  t6 = (-2.0 * (t3 * (-t40 + r) + t8 * ((t2 * (r - 0.1e1) + a52) * r - t39)) *
            r +
        t3 * (a52 * (-6.0) - 0.3e1 * epsi3) + (-t2 * t40 + t39) * epsi3 +
        4.0 * (epsi3 + t39) * t6) * SQR(t47) * SQR(t46);
  const Real& g_tt = t5;
  const Real g_pp = t43 * t17 * t14 * t48;
  const Real g_tp = -t4 * t9 * spin * t14 * t48;

  const Real gurr = t23 * t27 * t45 * t42;
  const Real guthth = r * t45;

  const Real dgttdr = t35 *
           (t28 * (t9 * (0.10e2 * r * t31 * t33 - 0.1e1) - t37 * r) +
            (-6.0) * (t2 * ((-0.5e1 / 0.3e1 + r) * r + t30) - t30 * t10 * t14) *
                t2 * t9);
  const Real dgtpdr = t38 * spin *
           ((0.20e2 * t9 * t31 * t33 + t29 * (-6.0)) * t4 / 0.2e1 -
            0.12e2 * r *
                ((t24 / 0.6e1 + t2 / 0.3e1) * a13 + t3 +
                 t25 * (0.5e1 / 0.12e2 * t2 + t7 / 0.4e1)) *
                t9);
  const Real& dgrrdr = t6;
  //dgththdr = -epsi3 * pow(t42, 0.2e1) + t22;
  const Real dgththdr = -epsi3 * SQR(t42) + t22;
  const Real dgppdr =
      0.10e2 * t38 * t9 *
          (-t31 * t17 * t33 * t42 + t20 * t32 -
           0.4e1 / 0.5e1 * t3 *
               ((-0.7e1 / 0.4e1 + r) * r + 0.3e1 / 0.4e1 * t7) * t7 * t14) +
      t38 * t17 * t42 * (-t43 + t37);

  const Real dgttdth = 2.0 * t18 * t1 * t11 * t35 * (r * t28 + t13 * t9) -
            4.0 * t2 * t1 * t7 * t10 * t28 * t11 * t9 * t36;
  const Real dgtpdth = -t19;
  const Real dgrrdth = -t21 * t2 * t11 * t47 * t46;
  const Real dgththdth = -t21 * t11;
  const Real& dgppdth = t44;

  const Real hgurr = 0.5 * gurr;
  const Real hguthth = 0.5 * guthth;

  const Real ch_rtt = -hgurr * dgttdr;
  const Real ch_rtp = -hgurr * dgtpdr;
  const Real ch_rrr = hgurr * dgrrdr;
  const Real ch_rrth = hgurr * dgrrdth;
  const Real ch_rthth = -hgurr * dgththdr;
  const Real ch_rpp = -hgurr * dgppdr;
  const Real ch_thtt = -hguthth * dgttdth;
  const Real ch_thtp = -hguthth * dgtpdth;
  const Real ch_thrr = -hguthth * dgrrdth;
  const Real ch_thrth = hguthth * dgththdr;
  const Real ch_ththth = hguthth * dgththdth;
  const Real ch_thpp = -hguthth * dgppdth;

  const Real denom = (g_tt * g_pp - g_tp * g_tp);

  const Real kt = -(g_pp + b * g_tp) / denom;
  const Real kphi = (g_tp + b * g_tt) / denom;

  diffs[0] = vars[3];
  diffs[1] = vars[4];
  diffs[2] = kphi;

  const Real kt2 = kt * kt;
  const Real kr2 = vars[3] * vars[3];
  const Real kth2 = vars[4] * vars[4];
  const Real kp2 = kphi * kphi;
  const Real ktp = kt * kphi;
  const Real krth = vars[3] * vars[4];

  diffs[3] = -(ch_rtt * kt2 + ch_rrr * kr2 + ch_rthth * kth2 + ch_rpp * kp2 +
               2.0 * (ch_rtp * ktp + ch_rrth * krth));
  diffs[4] = -(ch_thtt * kt2 + ch_thrr * kr2 + ch_ththth * kth2 +
               ch_thpp * kp2 + 2.0 * (ch_thtp * ktp + ch_thrth * krth));
}
