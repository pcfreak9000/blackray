#include "metric.hpp"

#include <cmath>

bool gluInvertMatrix(const Real m[16], Real invOut[16]) {
  Real inv[16], det;
  int i;

  inv[0] = m[5] * m[10] * m[15] - m[5] * m[11] * m[14] - m[9] * m[6] * m[15]
      + m[9] * m[7] * m[14] + m[13] * m[6] * m[11] - m[13] * m[7] * m[10];

  inv[4] = -m[4] * m[10] * m[15] + m[4] * m[11] * m[14] + m[8] * m[6] * m[15]
      - m[8] * m[7] * m[14] - m[12] * m[6] * m[11] + m[12] * m[7] * m[10];

  inv[8] = m[4] * m[9] * m[15] - m[4] * m[11] * m[13] - m[8] * m[5] * m[15]
      + m[8] * m[7] * m[13] + m[12] * m[5] * m[11] - m[12] * m[7] * m[9];

  inv[12] = -m[4] * m[9] * m[14] + m[4] * m[10] * m[13] + m[8] * m[5] * m[14]
      - m[8] * m[6] * m[13] - m[12] * m[5] * m[10] + m[12] * m[6] * m[9];

  inv[1] = -m[1] * m[10] * m[15] + m[1] * m[11] * m[14] + m[9] * m[2] * m[15]
      - m[9] * m[3] * m[14] - m[13] * m[2] * m[11] + m[13] * m[3] * m[10];

  inv[5] = m[0] * m[10] * m[15] - m[0] * m[11] * m[14] - m[8] * m[2] * m[15]
      + m[8] * m[3] * m[14] + m[12] * m[2] * m[11] - m[12] * m[3] * m[10];

  inv[9] = -m[0] * m[9] * m[15] + m[0] * m[11] * m[13] + m[8] * m[1] * m[15]
      - m[8] * m[3] * m[13] - m[12] * m[1] * m[11] + m[12] * m[3] * m[9];

  inv[13] = m[0] * m[9] * m[14] - m[0] * m[10] * m[13] - m[8] * m[1] * m[14]
      + m[8] * m[2] * m[13] + m[12] * m[1] * m[10] - m[12] * m[2] * m[9];

  inv[2] = m[1] * m[6] * m[15] - m[1] * m[7] * m[14] - m[5] * m[2] * m[15]
      + m[5] * m[3] * m[14] + m[13] * m[2] * m[7] - m[13] * m[3] * m[6];

  inv[6] = -m[0] * m[6] * m[15] + m[0] * m[7] * m[14] + m[4] * m[2] * m[15]
      - m[4] * m[3] * m[14] - m[12] * m[2] * m[7] + m[12] * m[3] * m[6];

  inv[10] = m[0] * m[5] * m[15] - m[0] * m[7] * m[13] - m[4] * m[1] * m[15]
      + m[4] * m[3] * m[13] + m[12] * m[1] * m[7] - m[12] * m[3] * m[5];

  inv[14] = -m[0] * m[5] * m[14] + m[0] * m[6] * m[13] + m[4] * m[1] * m[14]
      - m[4] * m[2] * m[13] - m[12] * m[1] * m[6] + m[12] * m[2] * m[5];

  inv[3] = -m[1] * m[6] * m[11] + m[1] * m[7] * m[10] + m[5] * m[2] * m[11]
      - m[5] * m[3] * m[10] - m[9] * m[2] * m[7] + m[9] * m[3] * m[6];

  inv[7] = m[0] * m[6] * m[11] - m[0] * m[7] * m[10] - m[4] * m[2] * m[11]
      + m[4] * m[3] * m[10] + m[8] * m[2] * m[7] - m[8] * m[3] * m[6];

  inv[11] = -m[0] * m[5] * m[11] + m[0] * m[7] * m[9] + m[4] * m[1] * m[11]
      - m[4] * m[3] * m[9] - m[8] * m[1] * m[7] + m[8] * m[3] * m[5];

  inv[15] = m[0] * m[5] * m[10] - m[0] * m[6] * m[9] - m[4] * m[1] * m[10]
      + m[4] * m[2] * m[9] + m[8] * m[1] * m[6] - m[8] * m[2] * m[5];

  det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

  if (det == 0)
    return false;

  det = 1.0 / det;

  for (i = 0; i < 16; i++)
    invOut[i] = inv[i] * det;

  return true;
}

void metric(Real r, Real th, Real mn[][4]) {
  Real t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15,
      t16;
  Real g_tt, g_rr, g_thth, g_pp, g_tp; /*metric components with lowered indicies*/

  t1 = cos(th);
  t2 = spin * spin;
  t3 = r * r;
  //t4 = pow(t3, 0.2e1);
  t4 = SQR(t3);
  t5 = t3 * t4;
  //t1 = t2 * pow(t1, 0.2e1);
  t1 = t2 * SQR(t1);
  t6 = (t1 + t3) * r + epsi3;
  t7 = a22 + t3;
  t8 = sin(th);
  t9 = t3 + t2;
  //t8 = pow(t8, 0.2e1);
  t8 = SQR(t8);
  t10 = r * t3 + a13;
  t11 = -t2 * r * t7 * t8 + t10 * t9;
  t12 = -0.2e1 * r + t9;
  t13 = a52 + t3;
  t14 = 0.1e1 / r;
  t15 = t2 * a22;
  t16 = 0.1e1 / t12;
  t13 = 0.1e1 / t13;
  //t11 = pow(t11, -0.2e1);
  t11 = 1.0 / SQR(t11);

  g_tt = t6 * (t2 * SQR(t7) * t8 + 2.0 * r * t4 - t4 * t9) * r * t11;
  g_rr = t6 * r * t16 * t13;
  g_thth = t14 * epsi3 + t1 + t3;
  g_pp = (SQR(t9) * SQR(t10) - t2 * t5 * t12 * t8) * t6 * t8 * t11
      * t14;
  g_tp = -(2.0 * t5
      + t3 * (a13 * (a22 + t2) + ((r * a22 + a13) * r + t15) * r) + t15 * a13)
      * spin * t6 * t8 * t11;

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      mn[i][j] = 0.0;
    }
  }

  mn[0][0] = g_tt;
  mn[0][3] = g_tp;
  mn[1][1] = g_rr;
  mn[2][2] = g_thth;
  mn[3][0] = mn[0][3];
  mn[3][3] = g_pp;
}
void metric_inv(Real r, Real th, Real mn[][4]) {
  metric(r, th, mn);
  Real mn_1d[16];
  Real mn_1d_inv[16];
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      mn_1d[i + j * 4] = mn[i][j];
    }
  }
  gluInvertMatrix(mn_1d, mn_1d_inv);
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      mn[i][j] = mn_1d_inv[i + j * 4];
    }
  }
}

void metric_rderivatives(Real r, Real th, Real dmn[][4]) {
  Real t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15,
      t16, t17, t18, t19, t20, t21, t22, t23;
  Real dgttdr, dgtpdr, dgppdr;

  t1 = cos(th);
  t2 = spin * spin;
  t3 = r * r;
  //t4 = pow(t3, 0.2e1);
  t4 = SQR(t3);
  t5 = t3 * t4;
  //t1 = t2 * pow(t1, 0.2e1);
  t1 = t2 * SQR(t1);
  t6 = (t1 + t3) * r + epsi3;
  t7 = a22 + t3;
  t8 = sin(th);
  t9 = t3 + t2;
  t10 = -0.2e1 * r + t9;
  //t8 = pow(t8, 0.2e1);
  t8 = SQR(t8);
  t11 = r * t3 + a13;
  t12 = t9 * t11;
  t13 = -t2 * r * t7 * t8 + t12;
  t1 = t1 / 0.3e1 + t3;
  t14 = 0.2e1 / 0.3e1 * t2;
  t15 = 0.3e1 / 0.5e1 * t2;
  t16 = (t15 + t3) * r + 0.2e1 / 0.5e1 * a13;
  t15 = r * t16 - t15 * (a22 / 0.3e1 + t3) * t8;
  t13 = 0.1e1 / t13;
  //t17 = pow(t13, 0.2e1);
  t17 = SQR(t13);
  t18 = 0.3e1 * t1;
  t19 = a22 + t2;
  t20 = r * a22;
  t21 = t2 * a22;
  t22 = 0.1e1 / 0.2e1;
  t23 = t17 * t8;
  //t9 = -t2 * t5 * t10 * t8 + pow(t9, 0.2e1) * pow(t11, 0.2e1);
  t9 = -t2 * t5 * t10 * t8 + SQR(t9) * SQR(t11);
  t11 = 0.1e1 / r;

  dgttdr =
      t17 * ((t6 * (0.10e2 * r * t15 * t13 - 0.1e1) - t18 * r) *
//                 (-t2 * pow(t7, 0.2e1) * t8 + t10 * t4) -
          (-t2 * SQR(t7) * t8 + t10 * t4) -
             0.6e1 * t6 *
                 (t3 * ((-0.5e1 / 0.3e1 + r) * r + t14) - t14 * t7 * t8) * t3);
  dgtpdr =
      t23 * spin *
      ((0.20e2 * t6 * t15 * t13 - 0.6e1 * t1) *
           (t22 * (t3 * (a13 * t19 + ((a13 + t20) * r + t21) * r) + t21 * a13) +
            t5) -
       0.12e2 * t6 *
           ((t19 / 0.6e1 + t3 / 0.3e1) * a13 + t4 +
            t20 * (0.5e1 / 0.12e2 * t3 + t2 / 0.4e1)) *
           r);
  dgppdr = 0.10e2 * t23 * t6 *
               (-t15 * t9 * t11 * t13 + t12 * t16 -
                0.4e1 / 0.5e1 * t4 *
                    ((-0.7e1 / 0.4e1 + r) * r + 0.3e1 / 0.4e1 * t2) * t2 * t8) +
           t23 * t9 * t11 * (-t11 * t6 + t18);
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      dmn[i][j] = 0.0;
    }
  }
  dmn[0][0] = dgttdr;
  dmn[0][3] = dgtpdr;
  dmn[3][0] = dmn[0][3];
  dmn[3][3] = dgppdr;
}
