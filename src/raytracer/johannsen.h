// #include "metric.h"

const  double Piby2 = 1.5707963267948966192;
double a, M;



void cache();
void rcache();

double R(double vars[5]);

double kerr_rms(double a);
double find_hitting_angle(double nvars[5], double ovars[5]);
double find_radii(double nvars[5], double ovars[5]);
double interp_lin_1d(double ifac_r, double rlo, double rhi);
double cal_q(double delta, double h);
double cal_delta_inc(double delta, double theta, double h, double r, const Real vars[], const Real prevVars[]);
double cal_kt(double r, double theta,double b);
double cal_gi(double delta, double theta, double h, double r);
double r0_grid(int i);
double r_mean(int i);
double I0(int i);
double II(double r);
double flat_I(double r,  double h);
void sortd(double ds[]);



// extern double a13;

void metric(double r, double th, double g[][4]);
void metric_rderivatives(double r, double th, double dg[][4]);
void metric_r2derivatives(double r, double th, double dg2[][4]);
void uppermetric(double r, double th, double rth[]);
double find_isco();

// void find_isco( double spin,  double spin2,  double eps3,  double a13,  double a22,  double a52,  double z1,  double& isco);



void metric(double r, double th, double g[][4])
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19;
  double gtt, grr, gthth, gpp, gtp;
  double spin = chi;

  t1 = r * r;
  t2 = pow(t1, 0.2e1);
  t3 = r * t1;
  t4 = t1 * t2;
  t5 = a22 + t1;
  t6 = sin(th);
  t6 = pow(t6, 0.2e1);
  t7 = spin * spin;
  t8 = cos(th);
  t8 = t7 * pow(t8, 0.2e1);
  t9 = (t8 + t1) * r + eps3;
  t10 = t7 + t1;
  t11 = a13 + t3;
  t12 = -t7 * r * t5 * t6 + t10 * t11;
  t13 = 0.1e1 / r;
  t8 = eps3 * t13 + t1 + t8;
  t14 = 0.2e1 * r;
  t15 = -t14 + t10;
  t16 = pow(t13, 0.2e1);
  t17 = a52 * t16 + 0.1e1;
  t18 = a22 * t16 + 0.1e1;
  t16 = t10 * (a13 * t13 * t16 + 0.1e1);
  t19 = -t7 * t18 * t6 + t16;
  t15 = 0.1e1 / t15;
  t17 = 0.1e1 / t17;
  t19 = pow(t19, -0.2e1);
  t12 = pow(t12, -0.2e1);

  gtt = (0.2e1 * r * t2 - t7 * (-t6 * pow(t5, 0.2e1) + t2) - t4) * r * t9 * t12;
  grr = t8 * t15 * t17;
  gthth = t8;
  gpp = t6 * (t6 * (-t4 * t7 * t10 + 0.2e1 * t3 * t2 * t7) + pow(t10, 0.2e1) * pow(t11, 0.2e1)) * t9 * t12 * t13;
  gtp = -spin * (t16 * t18 - t1 + t14 - t7) * t8 * t6 * t19;

  g[0][0] = gtt;
  g[0][3] = gtp;
  g[1][1] = grr;
  g[2][2] = gthth;
  g[3][0] = g[0][3];
  g[3][3] = gpp;
}

void metric_rderivatives(double r, double th, double dg[][4])
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29;
  double dgttdr, dgtpdr, dgppdr;
  double spin = chi;

  t1 = r * r;
  t2 = pow(t1, 0.2e1);
  t3 = r * t1;
  t4 = t1 * t2;
  t5 = a22 + t1;
  t6 = sin(th);
  t6 = pow(t6, 0.2e1);
  t7 = spin * spin;
  t8 = t5 * t6;
  t9 = cos(th);
  t9 = t7 * pow(t9, 0.2e1);
  t10 = (t1 + t9) * r + eps3;
  t11 = t1 + t7;
  t12 = a13 + t3;
  t13 = t11 * t12;
  t14 = t8 * t7 * r - t13;
  t15 = 2;
  t5 = (double) t15 * r * t2 - t7 * (-t6 * pow(t5, 0.2e1) + t2) - t4;
  t16 = 0.3e1 * t1;
  t17 = t16 + t9;
  t18 = (0.5e1 / 0.2e1 * t1 + 0.3e1 / 0.2e1 * t7) * r + a13;
  t16 = r * t18 - t7 * (a22 + t16) * t6 / 0.2e1;
  t14 = 0.1e1 / t14;
  t19 = pow(t14, 0.2e1);
  t4 = t6 * ((double) t15 * t3 * t2 * t7 - t4 * t7 * t11) + pow(t11, 0.2e1) * pow(t12, 0.2e1);
  t12 = 0.1e1 / r;
  t20 = t10 * t12;
  t21 = t7 * a22;
  t22 = 0.3e1 / 0.5e1;
  t23 = pow(t12, 0.2e1);
  t24 = t12 * t23;
  t25 = a13 * t24 + 0.1e1;
  t26 = a22 * t23 + 0.1e1;
  t27 = t11 * t25;
  t28 = (double) t15 * r;
  t29 = -t7 * t26 * t6 + t27;
  t29 = 0.1e1 / t29;

  dgttdr = 0.4e1 * t10 * r * t19 * (r * (t3 * (-0.3e1 / 0.2e1 * r + 0.5e1 / 0.2e1) - t7 * (t1 - t8)) + t16 * t5 * t14) + t5 * t19 * (r * t17 + t10);
  dgppdr = 0.4e1 * t6 * t10 * t19 * (t13 * t18 + (0.7e1 / 0.2e1 * (-0.4e1 / 0.7e1 * r + 0.1e1) * r - 0.3e1 / 0.2e1 * t7) * t7 * t2 * t6 + t4 * t16 * t14 * t12) + t6 * t4 * t19 * t12 * (t17 - t20);
  dgtpdr = t6 * spin * (0.5e1 * t20 * (a13 * (t1 * (t1 / 0.5e1 + (a22 + t7) * t22) + t21) - 0.2e1 / 0.5e1 * t3 * (t3 - t21)) * t19 + (t27 * t26 - t1 + t28 - t7) * pow(t29, 0.2e1) * (eps3 * t23 - t28 + (double) t15 * (eps3 * t12 + t1 + t9) * t29 * ((double) t15 * (t21 * t24 * t6 + r * t25) - 0.3e1 * t11 * a13 * pow(t23, 0.2e1))));

  dg[0][0] = dgttdr;
  dg[0][3] = dgtpdr;
  dg[3][0] = dg[0][3];
  dg[3][3] = dgppdr;
}

void metric_r2derivatives(double r, double th, double dg2[][4])
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42, t43, t44, t45, t46, t47, t48, t49, t50, t51, t52, t53, t54, t55, t56, t57, t58, t59, t60, t61, t62, t63, t64, t65, t66, t67, t68, t69, t70, t71, t72, t73, t74, t75, t76, t77, t78, t79, t80, t81, t82, t83, t84, t85, t86, t87, t88, t89, t90, t91, t92, t93, t94, t95, t96, t97, t98, t99, t100, t101, t102, t103, t104, t105, t106, t107, t108, t109, t110, t111;
  double dgttdr2, dgtpdr2, dgppdr2;
  double spin = chi;

  t1 = r * r;
  t2 = pow(t1, 0.2e1);
  t3 = pow(t2, 0.2e1);
  t4 = r * t2;
  t5 = r * t1;
  t6 = t5 * t3;
  t7 = r * t3;
  t8 = t5 * t2;
  t9 = t1 * t2;
  t10 = a22 + t1;
  t11 = eps3 + t5;
  t12 = sin(th);
  t13 = 2;
  t14 = (double) t13 * a22;
  t15 = spin * spin;
  t16 = pow(t15, 0.2e1);
  t17 = t15 * t16;
  t18 = 0.3e1 * t15;
  t19 = a22 * a22;
  t20 = a22 * t19;
  t21 = 0.5e1 * t15;
  t22 = 0.6e1 * a22;
  t23 = a22 - t18;
  t24 = t5 * a22;
  t25 = a13 + t5;
  t26 = 0.3e1 / 0.2e1;
  t27 = 0.1e1 / 0.2e1;
  t28 = t19 * t26;
  t29 = a13 * a22;
  t30 = 0.4e1 * a13;
  t31 = cos(th);
  t32 = t26 * t15;
  t33 = t22 - t32;
  t34 = a13 + eps3;
  t35 = 0.4e1 * t15;
  t36 = 0.12e2 * a22;
  t37 = a13 / 0.7e1;
  t38 = 0.6e1 / 0.7e1;
  t39 = 0.3e1 * a22;
  t40 = 0.6e1 * eps3;
  t41 = a13 * eps3;
  t42 = (double) t13 * a13;
  t43 = eps3 + t42;
  t44 = a22 * eps3;
  t45 = -0.21e2 / 0.2e1 * eps3;
  t46 = t29 * eps3;
  t47 = t15 * t20;
  t31 = pow(t31, 0.2e1);
  t48 = t15 * t31;
  t49 = (double) t13 * t15;
  t50 = 36;
  t51 = -0.26e2 / 0.3e1;
  t52 = 0.2e1 / 0.3e1;
  t53 = t52 * a22;
  t54 = 0.12e2 * a13;
  t55 = 0.13e2 / 0.3e1;
  t56 = 0.5e1 / 0.2e1 * t19;
  t57 = t29 * t15;
  t58 = 0.4e1 * a22;
  t59 = a13 * a13;
  t60 = a13 * t59;
  t61 = t59 * t16;
  t62 = t61 * t19;
  t63 = 0.3e1 * t59;
  t64 = 0.16e2 * a13;
  t65 = t59 * a22;
  t66 = t26 * eps3;
  t67 = 0.10e2 / 0.3e1;
  t68 = 0.7e1 / 0.2e1;
  t69 = t27 * t15;
  t70 = t67 * a13;
  t71 = t52 * eps3;
  t72 = 0.3e1 * eps3;
  t73 = t27 * t19;
  t74 = 0.28e2 / 0.3e1 * a22;
  t75 = (double) t13 * eps3;
  t76 = t27 * eps3;
  t77 = -a13 * t52 - t76;
  t78 = 0.15e2 / 0.2e1;
  t79 = 0.8e1 / 0.3e1;
  t80 = t79 * t15;
  t81 = eps3 / 0.4e1;
  t82 = a13 + t81;
  t83 = t19 * eps3;
  t84 = a13 - 0.8e1 / 0.7e1 * eps3;
  t85 = 0.1e1 / 0.3e1;
  t86 = t85 * t15;
  t87 = 0.18e2 * a22;
  t88 = a13 + t75;
  t89 = a22 * t88;
  t90 = a22 - t15;
  t91 = t59 * eps3;
  t92 = t52 * t15;
  t93 = t3 * ((r * t27 - t79) * r + t92) - t30 * (0.7e1 / 0.6e1 * t16 * t84 + t83 + t86 * a22 * (a13 - 0.13e2 / 0.2e1 * eps3));
  t32 = t48 * (t1 * (((t1 * (t5 * ((double) t13 * t4 + 0.6e1 * t15 * (a13 * t51 + t53) + 0.6e1 * t29) + t35 * (a13 * (a13 - 0.14e2 * a22) + t56 * t15)) + t63 * (t15 * (-t58 + t21) + t19)) * r - t64 * t19 * t16) * r - 0.8e1 * t65 * (a22 - t32) * t15) + 0.12e2 * (a13 * (t15 * (-a22 * t55 - t49) + t19) + t2 * (-a13 + 0.16e2 / 0.3e1 * a22 + t1)) * t8 + t4 * (t5 * (t1 * (t1 * (t36 - t49) + a22 * (a22 * (double) t50 - t35)) + a13 * (a13 - 0.20e2 * a22) + t15 * (0.34e2 * t19 + t54)) + t57 * (-0.28e2 * a22 - 0.18e2 * t15)) + t62) / 0.6e1;
  t50 = (int) (0.2e1 / 0.7e1);
  t94 = 0.3e1 / 0.7e1;
  t95 = (double) t50 * t15;
  t96 = 0.1e1 / 0.42e2;
  t97 = -0.23e2 / 0.42e2;
  t98 = t15 * a13;
  t99 = 0.4e1 / 0.3e1 * t15;
  t100 = 0.1e1 / 0.5e1;
  t101 = eps3 / 0.14e2 + a13;
  t102 = 0.4e1 / 0.5e1;
  t103 = a13 * (a13 - 0.8e1 / 0.5e1 * eps3);
  t104 = a13 - 0.26e2 / 0.17e2 * eps3;
  t105 = t85 * eps3;
  t106 = a13 - t105;
  t107 = 0.4e1 / 0.15e2;
  t108 = a13 * t1;
  t31 = r * (0.7e1 / 0.5e1 * t31 * (t1 * (t2 * ((t5 * ((-r / 0.14e2 + (double) t50) * r - t15 / 0.21e2) + 0.4e1 / 0.3e1 * t98) * r - t37 * (a13 + 0.16e2 * t15)) - 0.16e2 / 0.21e2 * t61) + (t1 * (((t5 * ((a13 * t94 - t95) * r - a13 * t38 + t16 * t96) + a13 * (0.25e2 / 0.21e2 * t16 + t37)) * r + t98 * (a13 * t97 - t95)) * r + t98 * (t16 * (double) t50 + 0.4e1 / 0.7e1 * a13)) + t61) * r - 0.5e1 / 0.14e2 * t17 * t59) + t108 * ((t16 * (-a13 * t68 + t75) + t41) * t107 - 0.9e1 / 0.5e1 * t2 * (-0.4e1 / 0.9e1 * t15 + a13 - 0.10e2 / 0.9e1 * eps3))) + t15 * (-0.19e2 / 0.30e2 * t91 * t1 - t91 * t86);
  t50 = 60;
  t12 = pow(t12, 0.2e1);
  t97 = pow(t12, 0.2e1);
  t107 = pow(t10, 0.2e1);
  t109 = (double) t13 * r;
  t20 = (-0.4e1 * (double) t13 * t1 * ((t3 * (t34 - t22) + t47 * t34) * r + t41 * t19 * t23) + 0.12e2 * t4 * t3 - 0.28e2 * t2 * (t4 * (a22 * (eps3 * t38 + a13 - t39) + t15 * (0.3e1 / 0.14e2 * eps3 + t37)) + t46 * (a22 - 0.13e2 / 0.7e1 * t15)) + 0.4e1 * (t2 * (t19 * (a22 * t34 + t15 * (-0.9e1 * a13 + t45)) + t5 * (t1 * (t1 * (t33 + t1) + a22 * (t36 - t35)) + ((0.13e2 * a22 - 0.9e1 / 0.2e1 * t15) * a22 - t40) * a22 + t41)) + t48 * (t1 * (t1 * (t5 * (t2 * t27 + t28) - 0.3e1 * t29 * t23) + t30 * t19 * t15) + (double) t13 * t9 * (-a22 * (a13 - 0.5e1 * a22) + t18 * a13 + t24) + t4 * (t5 * (a13 - t14) + t19 * (t22 - t21)) + t25 * t15 * t20)) * r + 0.4e1 * t47 * t41 - 0.16e2 * t9 * (t15 * (-t42 * eps3 - t20) + t44 * (-a22 * t26 + a13) + (a22 * (a13 + 0.9e1 / 0.4e1 * eps3) + t15 * t43) * a22 * r)) * t16 * t97;
  t23 = (0.12e2 * (t4 * t93 - t85 * t2 * (t5 * (t15 * (0.11e2 * a13 * (a13 - 0.54e2 / 0.11e2 * a22 - 0.30e2 / 0.11e2 * eps3) + 0.6e1 * t44 + t28 * t15) + t29 * (a13 + 0.7e1 * eps3)) + t41 * (-t15 * (a13 + t87) + t29) - 0.9e1 * t19 * (a13 - 0.5e1 / 0.6e1 * eps3) * t16) - 0.4e1 / 0.3e1 * t9 * (t2 * (a22 * (a13 + t72) + t15 * (-a13 * t26 + eps3 - t39)) - t15 * (0.6e1 * (-t76 + t19) * a13 - 0.6e1 * t83 - t14 * t82 * t15) - t46) + t1 * (t9 * ((t1 * ((a22 * t51 + t49 - t66) * r + t15 * (t14 - t69) - t19 * t68 + t70 - t71) - t15 * (t53 * t15 - 0.4e1 * a13 + 0.4e1 * t73) - a13 * (a13 - t75 - t74) - 0.32e2 / 0.3e1 * t44) * r + (-t44 * t78 + t80 * (a13 - 0.3e1 / 0.4e1 * eps3)) * a22 + (-t71 + t19) * a13 + t16 * t77) - t91 * (a22 - 0.5e1 / 0.3e1 * t15) * t90) + t57 * (t5 * (t15 * (eps3 * t55 - 0.5e1 * a13) + t89) + t41 * t90) - a13 * t16 * t19 * (a13 - t75) * r) * r - 0.12e2 * t32) * t15 * t12;
  t28 = t15 + t1;
  t32 = t15 * r;
  t10 = t32 * t10 * t12 - t25 * t28;
  t34 = 0.3e1 / 0.10e2 * t15;
  t36 = -0.7e1 / 0.3e1;
  t37 = t15 * a22;
  t38 = 0.20e2 / 0.3e1 * a22;
  t46 = 0.3e1 / 0.4e1 * t15;
  t47 = 0.21e2;
  t51 = 0.14e2 / 0.3e1;
  t53 = -0.1e1 / 0.8e1;
  t55 = 0.13e2 / 0.6e1;
  t93 = eps3 / 0.16e2;
  t105 = a13 + t105;
  t110 = t35 * a22;
  t72 = -a13 + t72;
  t111 = t68 * t15;
  t33 = t48 * t27 * (t2 * (t1 * (t5 * (t1 * (-t5 * t52 + 0.4e1 / 0.3e1 * a13 + 0.64e2 / 0.3e1 * a22) + (t95 * (a22 + a13) + t29) * t51) - t79 * (t15 * (-t46 * t19 - 0.5e1 / 0.2e1 * a13 * (a13 - 0.14e2 / 0.5e1 * a22)) + t65)) + t59 * (t15 * (t15 * t47 - t58) + t19) * t85) + 0.4e1 * t1 * (t61 * a22 + t6) + t4 * (t1 * ((t1 * (-0.4e1 / 0.9e1 * t1 * t33 + a22 * (-0.20e2 / 0.3e1 * t15 + t74)) + a13 * (a13 - t38) - 0.4e1 * t15 * (-a13 - t56 + t37)) * r + t29 * (t35 + t38)) + t57 * (t58 + t49)) + t62);
  t33 = t5 * (t9 * (((t1 * (0.4e1 * a22 * t55 + 0.4e1 * eps3 * t53 + t1 * (r * t26 + t79) + 0.4e1 * a13 - 0.4e1 * t69) + 0.32e2 / 0.3e1 * a22 * (a13 - 0.5e1 / 0.16e2 * eps3) + 0.32e2 / 0.3e1 * t15 * (0.7e1 / 0.16e2 * a13 - 0.3e1 / 0.8e1 * a22 + t93)) * r - 0.28e2 / 0.3e1 * a22 * t84 - 0.28e2 / 0.3e1 * t15 * (t94 * (t19 * t36 + a13) - t37 / 0.7e1)) * r + 0.7e1 * (t44 * t27 + 0.40e2 / 0.21e2 * t15 * (a13 - 0.9e1 / 0.20e2 * eps3)) * a22 + 0.7e1 * (0.2e1 / 0.21e2 * eps3 + t19) * a13 - t16 * t77) + t99 * t91 * (a22 + t111)) + t33;
  t38 = -0.77e2 / 0.6e1;
  t53 = t65 * t17;
  t55 = 0.15e2 / 0.4e1;
  t58 = 0.3e1 * a13;
  t74 = 0.17e2 / 0.4e1;
  t77 = t15 * t59;
  t35 = t3 * (a13 * (t17 * (double) t13 * t43 + t91 + t98 * (t37 * t47 - t75)) + (t4 * (t1 * (((-a22 * t27 + t15 * t55 + t1) * r + t35 + t58 - t76) * r + t15 * (t14 + t111) + 0.8e1 * a13 - 0.5e1 * eps3) + t15 * (t15 * (a22 * t68 + t46) + 0.14e2 * t101) + (0.33e2 / 0.4e1 * a13 - 0.3e1 * eps3) * a13) + t15 * (t15 * (-0.18e2 * a13 * t106 + 0.6e1 * t37 * t106) - t63 * (a13 - 0.8e1 * eps3)) - t65 * t72) * r) + t44 * t17 * t60;
  t6 = t48 * t26 * (t1 * (((-t27 * t6 - a13 * (t16 * (a13 * t47 - t18 * a22) + t59 * t90) * t85) * r + 0.7e1 * t61 * (a22 + t69)) * r + t92 * t60 * t90) - (double) t13 * (t6 * (a13 - t15 + t1) - t53) * r + 0.4e1 * a13 * t4 * (-t4 * (a22 + t86) + t98 * (t15 * t67 + a22)) + t9 * (((-a13 * (t15 * (t39 - t92) + a13) + (t1 * (t1 * (-a22 * t67 - t86) + t15 * (-0.19e2 / 0.3e1 * a22 + t15 / 0.6e1) + 0.6e1 * a13) + t15 * (-t110 + t64) + t63) * r) * r - a22 * (t59 + t17) - t98 * (a13 * t38 - t49)) * r + t98 * (t14 * t15 - t30)) + t16 * t60 * t90);
  t6 = t7 * ((-0.5e1 * t16 * (-0.11e2 / 0.5e1 * a13 - 0.2e1 / 0.5e1 * eps3) - 0.5e1 * t103 + 0.21e2 * t37 * (a13 - 0.13e2 / 0.21e2 * eps3) + (0.9e1 / 0.2e1 * a22 * (a13 - 0.13e2 / 0.9e1 * eps3) + 0.9e1 / 0.2e1 * t15 * (t79 * (a13 + t93) - t86)) * t1) * r + 0.6e1 * t15 * (t59 * t74 + a22 * t16 / 0.6e1) - 0.6e1 * a13 * (-a22 * (a13 - t66) + t16)) + t6;
  t18 = pow(t28, 0.2e1);
  t25 = pow(t25, 0.2e1);
  t14 = (0.12e2 * ((t1 * (t68 * (a22 * (a22 + eps3 / 0.21e2) - t96 * t15 * eps3) + t14 * t1) - 0.5e1 / 0.6e1 * t44 * (-t102 * t15 + a22)) * r + 0.5e1 / 0.3e1 * t48 * (-t100 * t24 * (r + 0.1e1) - t9 / 0.20e2 + r * (r * a22 * (-0.3e1 / 0.4e1 * a22 + t34) + t19) - t34 * t19)) * r + 0.12e2 * t27 * t3 * (-r + 0.1e1) + 0.12e2 * (t1 * (t1 * (t1 * (-t15 / 0.12e2 - 0.11e2 / 0.6e1 * a22) + a22 * (a22 * t36 - t86)) + a22 * (-eps3 - 0.5e1 / 0.4e1 * t37)) + t83) * r - 0.12e2 * t81 * t15 * t19) * t15 * t9 * t12;
  t24 = t32 * t12;
  t14 = t24 * (0.6e1 * r * t33 + 0.6e1 * t85 * t57 * t1 * (t5 * (t15 * t72 + t89) + t41 * (a22 + 0.11e2 * t15)) + 0.6e1 * t9 * ((a13 * (t16 * (0.5e1 / 0.3e1 * a13 + t71) + t19 * t88 - t80 * t82 * a22) + ((t2 * (t1 * (a22 * t51 + t80) + t15 * (t22 + t69) + t19 * t78 - t70 + t71) + t15 * (a13 * (a13 + t71 - t87) + t44 * (double) t13) - t29 * t105 + t56 * t16) * r + t15 * (0.20e2 / 0.3e1 * (0.3e1 / 0.5e1 * eps3 + t19) * a13 + 0.20e2 / 0.3e1 * t73 * eps3 + t110 * (a13 - t71)) + t41 * (-0.4e1 / 0.3e1 * a22 + a13)) * r) * r - a22 * (-a22 * t27 * t43 * t16 + t91) + t21 * t41 * (a13 - 0.6e1 / 0.5e1 * a22)) + 0.6e1 * t62 * eps3 + t14);
  t19 = (t48 + t1) * r + eps3;
  t21 = a22 + t15;
  t22 = t5 * (-t37 + t5);
  t26 = a13 * (t1 * (t1 * t100 + 0.3e1 / 0.5e1 * t21) + t37) - 0.2e1 / 0.5e1 * t22;
  t33 = t37 * t12;
  t34 = 0.1e1 / r;
  t36 = pow(t34, 0.2e1);
  t38 = t34 * t36;
  t39 = a13 * t38 + 0.1e1;
  t44 = a22 * t36 + 0.1e1;
  t46 = t28 * t39;
  t47 = t46 * t44 - t1 + t109 - t15;
  t44 = -t15 * t44 * t12 + t46;
  t39 = (double) t13 * (r * t39 + t33 * t38) - t58 * t28 * pow(t36, 0.2e1);
  t46 = eps3 * t34 + t1 + t48;
  t10 = 0.1e1 / t10;
  t51 = pow(t10, 0.2e1);
  t55 = pow(t51, 0.2e1);
  t10 = t10 * t51;
  t44 = 0.1e1 / t44;
  t56 = pow(t44, 0.2e1);
  t57 = t12 * spin;
  t10 = -(double) t13 * t57 * (-((0.5e1 * t37 + t2) * a13 + 0.3e1 * t108 * t21 - (double) t13 * t22) * t19 * ((0.5e1 * t5 + t42) * r + 0.3e1 * t15 * t1 * (0.1e1 - t12) - t33) * t10 * t34 + t47 * t56 * (eps3 * t38 - t46 * t44 * (t38 * (t34 * (t54 * t28 * t34 - 0.6e1 * t33) - 0.10e2 * a13) + (double) t13) + 0.1e1)) + 0.5e1 * t51 * t34 * t26 * t12 * spin * (-t34 * (-(double) t13 * t5 + eps3 + t19) + 0.3e1 * t1 + t48) + 0.6e1 * t57 * (-t47 * t46 * pow(t56, 0.2e1) * pow(t39, 0.2e1) + (-(double) t13 * t2 + a13 * (t1 * t52 + t15) + (a13 + t32) * a22) * t19 * t51) + t57 * (-0.30e2 * t26 * ((t1 * t85 + t15) * a13 - t52 * r * (t33 + t2)) * t19 * t36 * t10 + 0.4e1 * t47 * (-eps3 * t36 + t109) * t44 * t56 * t39);
  t1 = (double) t13 * t12 * (pow(t18, 0.2e1) * pow(t25, 0.2e1) * t11 + t24 * (-0.4e1 * (double) t13 * t77 * t2 * (t41 * (a22 + t49) + t32 * (a22 * (a13 + t40) + t15 * (eps3 * t74 - t42))) + 0.4e1 * t27 * t5 * (t5 * (a13 * (-0.9e1 * t15 * (t105 * a22 * t16 + t91) + t41 * (0.15e2 * t16 + t29)) + t2 * (a13 * (t16 * (-0.85e2 / 0.2e1 * a13 - 0.10e2 * eps3) + t41 - 0.45e2 * t37 * (a13 - 0.7e1 / 0.15e2 * eps3)) + (t3 + t15 * (0.34e2 * a13 * t104 - 0.45e2 * t37 * (a13 - 0.17e2 / 0.45e2 * eps3)) + t17 * (-t66 - t30) + t59 * (a13 + t45)) * r)) + t53 * (a13 - 0.9e1 * eps3)) - 0.4e1 * t5 * t6 - 0.4e1 * t68 * t77 * t1 * (t4 * (a22 * (a13 + 0.9e1 / 0.7e1 * eps3) + t15 * (-0.13e2 / 0.7e1 * a13 + 0.109e3 / 0.14e2 * eps3)) + t41 * t15 * (a22 + t15 / 0.7e1)) - 0.4e1 * t35 + t14)) * t55 * t38;

  dgttdr2 = (t109 * t17 * pow(t107, 0.2e1) * t11 * t12 * t97 + t20 - t23 + (double) t50 * (t15 * t31 + t9 * (t15 * t8 - t17 * t85 * t43 - t91 + 0.68e2 / 0.3e1 * t98 * t104) / 0.10e2 - 0.2e1 / 0.5e1 * t2 * (-t3 * (a13 - t99 - t76) + t98 * (-0.6e1 * t15 * t106 + t41)) + r * (t3 * (((t16 / 0.15e2 - 0.16e2 / 0.15e2 * a13 + t71) * r + t15 * (t15 * t100 + a13 - 0.3e1 / 0.10e2 * eps3)) * r + t15 * (-t16 / 0.30e2 - 0.28e2 / 0.15e2 * t101) - t27 * a13 * (a13 - 0.6e1 / 0.5e1 * eps3)) + t61 * eps3) + t4 * (a13 * (t16 * (-0.67e2 / 0.2e1 * a13 + 0.29e2 * eps3) + t41) + t7) / 0.15e2 + t52 * (t16 * (a13 * t102 - eps3 * t100) + t103) * t3) * t5) * t55;
  dgppdr2 = t1;
  dgtpdr2 = t10;

  dg2[0][0] = dgttdr2;
  dg2[0][3] = dgtpdr2;
  dg2[3][0] = dg2[0][3];
  dg2[3][3] = dgppdr2;
}


void uppermetric(double r, double th, double rth[])
{
  double gurr, guthth;
  double t1, t2, t3, t4;
  double spin = chi;

  t1 = cos(th);
  t2 = 0.1e1 / r;
  t3 = r * r;
  t4 = spin * spin;
  t1 = t4 * pow(t1, 0.2e1) + eps3 * t2 + t3;
  t1 = 0.1e1 / t1;

  gurr = t1 * (-0.2e1 * r + t3 + t4) * (a52 * pow(t2, 0.2e1) + 0.1e1);
  guthth = t1;

  rth[0] = gurr;
  rth[1] = guthth;
}


/* Calculate ISCO radius */
double find_isco()
{
    double d2Veff, d2Veff2, E_var, Lz_var, Omega_var, denom, den, den_r, den_rr, num, num_r, num_rr;
    double m[4][4],dmdr[4][4],dmdr2[4][4];
    double l, j, r, d2Veff_last = 0, d2Veff_last2 = 1000, rin;
    double jmin, jmax, lmin, lmax, rmin, factor, rstep;

    rstep = 0.1;
    rmin = 1.1;
    jmax = 15.0;
    jmin = rmin;
    factor = 1.0e2;

    for(j=jmax; j>jmin && j>rmin; j-=rstep/factor)
    {
        r = j;

    /* Calculate 2nd-derivative of the effective potential - d2Veff */
        metric(r, Pi/2., m);
        metric_rderivatives(r, Pi/2., dmdr);
        metric_r2derivatives(r, Pi/2., dmdr2);

        Omega_var = (-dmdr[0][3] + sqrt(dmdr[0][3]*dmdr[0][3] - dmdr[0][0]*dmdr[3][3])) / dmdr[3][3];
        denom = sqrt(-(m[0][0] + 2.0*m[0][3]*Omega_var + m[3][3]*Omega_var*Omega_var));
        E_var = -(m[0][0] + m[0][3]*Omega_var) / denom;
        Lz_var = (m[0][3] + m[3][3]*Omega_var) / denom;

        den = m[0][3]*m[0][3]-m[0][0]*m[3][3];
        den_r = 2.0*m[0][3]*dmdr[0][3]-dmdr[0][0]*m[3][3]-m[0][0]*dmdr[3][3];
        den_rr = 2.0*dmdr[0][3]*dmdr[0][3]+2.0*m[0][3]*dmdr2[0][3]-dmdr2[0][0]*m[3][3]-2.0*dmdr[0][0]*dmdr[3][3]-m[0][0]*dmdr2[3][3];
        num = E_var*E_var*m[3][3]+2.0*E_var*Lz_var*m[0][3]+Lz_var*Lz_var*m[0][0];
        num_r = E_var*E_var*dmdr[3][3]+2.0*E_var*Lz_var*dmdr[0][3]+Lz_var*Lz_var*dmdr[0][0];
        num_rr = E_var*E_var*dmdr2[3][3]+2.0*E_var*Lz_var*dmdr2[0][3]+Lz_var*Lz_var*dmdr2[0][0];

        d2Veff = (num_rr + (-2.0*num_r*den_r - num*den_rr + 2.0*num*den_r*den_r/den) / den) / den;

    /* Search for when d2Veff flips sign */
        if( (d2Veff < 0.0 && d2Veff_last > 0.0) || (d2Veff > 0.0 && d2Veff_last < 0.0) )
        {
            rin = j;
            lmax = rin + rstep/factor;
            lmin = rin - rstep/factor;
      jmin = lmin;
            factor *= 100.0;
            d2Veff_last2 = 1000;

      /* Zoom in on where d2Veff flips sign to increase ISCO accuracy */
            while(factor < 1.0e10)
            {
                for(l=lmax; l>lmin && l>rmin; l-=rstep/factor)
                {
                    r = l;

          /* Calculate 2nd-derivative of the effective potential - d2Veff */
                    metric(r, Pi/2., m);
                    metric_rderivatives(r, Pi/2., dmdr);
                    metric_r2derivatives(r, Pi/2., dmdr2);

                    Omega_var = (-dmdr[0][3] + sqrt(dmdr[0][3]*dmdr[0][3] - dmdr[0][0]*dmdr[3][3])) / dmdr[3][3];
                    denom = sqrt(-(m[0][0] + 2.0*m[0][3]*Omega_var + m[3][3]*Omega_var*Omega_var));
                    E_var = -(m[0][0] + m[0][3]*Omega_var) / denom;
                    Lz_var = (m[0][3] + m[3][3]*Omega_var) / denom;

                    den = m[0][3]*m[0][3]-m[0][0]*m[3][3];
                    den_r = 2.0*m[0][3]*dmdr[0][3]-dmdr[0][0]*m[3][3]-m[0][0]*dmdr[3][3];
                    den_rr = 2.0*dmdr[0][3]*dmdr[0][3]+2.0*m[0][3]*dmdr2[0][3]-dmdr2[0][0]*m[3][3]-2.0*dmdr[0][0]*dmdr[3][3]-m[0][0]*dmdr2[3][3];
                    num = E_var*E_var*m[3][3]+2.0*E_var*Lz_var*m[0][3]+Lz_var*Lz_var*m[0][0];
                    num_r = E_var*E_var*dmdr[3][3]+2.0*E_var*Lz_var*dmdr[0][3]+Lz_var*Lz_var*dmdr[0][0];
                    num_rr = E_var*E_var*dmdr2[3][3]+2.0*E_var*Lz_var*dmdr2[0][3]+Lz_var*Lz_var*dmdr2[0][0];

                    d2Veff2 = (num_rr + (-2.0*num_r*den_r - num*den_rr + 2.0*num*den_r*den_r/den) / den) / den;

                    if(fabs(d2Veff2)<fabs(d2Veff_last2))
                        rin = l;

                    d2Veff_last2 = d2Veff2;
                }

                lmax = rin + rstep/factor;
                lmin = rin - rstep/factor;
                factor *= 100.0;
                d2Veff_last2 = 1000;
            }

            factor = 1.0e2;
        }

        d2Veff_last = d2Veff;
    }

    // isco = rin;
    return rin;
}





double kerr_rms(double a){
  //   accounts for negative spin
  double sign = 1.0;
  if (a<0) {
     sign = -1.0;
  }

  double Z1 = 1.0+pow(1.0-a*a,1.0/3.0)*(pow(1.0+a,1.0/3.0)+pow(1.0-a,1.0/3.0));
  double Z2=sqrt((3.0*a*a)+(Z1*Z1));

  return 3.0+Z2-sign*sqrt((3.0-Z1)*(3.0+Z1+(2*Z2)));
}



double flat_I(double r,  double h) {
  return h / pow(r*r + h*h, 3./2.) / (2 * M_PI);
  // return 1. / pow(pow(r/h, 2.) + 1., 3./2.) / (2. * M_PI * pow(h, 2.));
}


void sortd(double ds[]){
  int i, j;
  double tmp;
  for (i=0;i<3;i++) 
    for (j=i+1; j<3; j++)
      if (ds[i] > ds[j]) {
        tmp = ds[i];
        ds[i] = ds[j];
        ds[j] = tmp;
      }
  
}


double cal_q(double delta, double height) {
  double a = chi;
  // double q = sqrt(a*a - 2.0*h + h*h) * sin(delta) / sqrt(a*a + h*h);
  double q = (a*a + height*height) * sin(delta);
  q *= q;
  q /= (height*height - 2*height + a*a);
  q -= a*a;
  q = sqrt(fabs(q));
  // q = sqrt((h*h - 2*h + a*a) / (a*a + h*h)) * sin(delta);
  return q;
}

void redshift(double r, double th, double ktkp, double& gg);

double cal_delta_inc(double delta, double th, double height, double r, const Real vars[], const Real prevVars[]) {

  // th = M_PI/2.0;
  th = fabs(th);
  double theta = th;
  double rder,rder1,rder2,thder,thder1,thder2,beta,t1,c1,c2,g[4][4],dg[4][4], rth[2];
  double angle;
  double gg;
  double Zr, Zth, K, cs1, ss1, ss3, r3, gurr, guthth;
  
  metric(height, 0.0, g);
  Real g_tt_source = g[0][0];

  metric(r, th, g);
  double dg_dr[4][4];
  metric_rderivatives(r, th, dg_dr);
  uppermetric(r, th, rth);

  double Omega = (-dg_dr[0][3] + std::sqrt(dg_dr[0][3]*dg_dr[0][3] - dg_dr[0][0]*dg_dr[3][3])) / dg_dr[3][3];

  gg = sqrt(g_tt_source / (g[0][0] + 2 * g[0][3] * Omega + g[3][3] * Omega * Omega));
  if (std::isnan(gg)) {
    printf(" gg is %f, %f, %f, %f, %f\n", gg, -g_tt_source, (-g[0][0] - 2 * g[0][3] * Omega - g[3][3] * Omega * Omega), -g_tt_source / (-g[0][0] - 2 * g[0][3] * Omega - g[3][3] * Omega * Omega), Omega);
  }
  K = 3. / eta * Mdot;
  cs1 = cos(th);
  ss1 = sin(th);
  ss3 = ss1*ss1*ss1;
  r3 = r*r*r;
  // Zr = 0.5 * K * sqrt(johannsen_isco / (r3 * ss1)) - cs1;
  Zr = cs1 - 0.5 * K * sqrt(johannsen_isco / (r3 * ss1));
  Zth = - 0.5 * K * cs1 * sqrt(johannsen_isco / (r * ss3)) - r*ss1;

  gurr = rth[0];
  guthth = rth[1];

  double vr = (vars[3] + prevVars[3]) / 2.0;
  double vth = (vars[4] + prevVars[4]) / 2.0;
  vr = prevVars[3];
  vth = prevVars[4];

  // printf(" vr=%f, vth=%f\n", vr, vth);

  angle = 1. / sqrt(gurr * Zr * Zr + guthth * Zth * Zth) * (- gurr * g[1][1] * Zr * vr - guthth * g[2][2] * Zth * vth);
  angle *= sqrt(-(g[0][0] + 2 * g[0][3] * Omega + g[3][3] * Omega * Omega));
  angle = acos(angle);

  // printf(" %f %f \n", angle, acos(cal_q(delta, height) / (r * sqrt(-1 / (g[0][0] + 2 * g[0][3] * Omega + g[3][3] * Omega * Omega)))));
  

  if (std::isnan(angle))
    printf(" angle, %f, %f %f\n", gg, gurr * Zr * Zr + guthth * Zth * Zth, Zr * vr + Zth * vth);

  return angle;  
}

// double cal_gi(double delta, double theta, double h, double r) {
//   double a=chi;
//   return cal_kt(r, theta)*sqrt(h*h-2*h+a*a)/sqrt(h*h + a*a);
// }

double cal_kt(double r, double theta, double b) {
  double g_tt, g_pp, g_tp, denom, kt, a = chi, M=1.0;
  double r2, a2, s2, c2;
  r2 = r*r;
  a2 = a*a;
  s2 = sin(theta)*sin(theta);
  c2 = cos(theta)*cos(theta);
  g_tt = -1.0 + (2.0*M*r)/(r2 + a2*c2);
  // g_tt = -1. + (2.*M*r)/(pow(r,2) + pow(a,2)*pow(cos(theta),2));
  // g_pp = pow(sin(theta),2)*(pow(a,2) + pow(r,2) + (2*pow(a,2)*M*r*pow(sin(theta),2))/(pow(r,2) + pow(a,2)*pow(cos(vars[1]),2)));
  g_pp = s2*(a2 + r2 + (2.0*a2*M*r*s2)/(r2 + a2*c2));
  g_tp = (-2.0*a*M*r*s2)/(r2 + a2*c2);
  // g_tp = (-2.0*a*M*r*pow(sin(theta),2))/(pow(r,2) + pow(a,2)*pow(cos(theta),2));  

  
  denom = (g_tt*g_pp-g_tp*g_tp);
  kt = -(g_pp+b*g_tp)/denom;
  // printf(" gtt=%f, gpp=%f, gtp=%f ", g_tt, g_pp, g_tp);
  return kt;
}


//void cache() {
//  int i = 0;
//  for (;i<N;i++)
//    prevVars[i] = vars[i];
//}
//
//void rcache() {
//  int i = 0;
//  for (;i<N;i++)
//    vars[i] = prevVars[i];
//}

double R(double vars[5]) {
  double x, y, z;
  
  x = vars[0]*sin(vars[1])*cos(vars[4]);
  y = vars[0]*sin(vars[1])*sin(vars[4]);
  z = vars[0]*cos(vars[1]);

  return sqrt(x*x + y*y + z*z);
}

double interp_lin_1d(double ifac_r, double rlo, double rhi){
  return ifac_r*rhi + (1.0-ifac_r)*rlo;
}


double find_radii(double nvars[5], double ovars[5]) {
  return (ovars[0]+nvars[0])/2.0;//sqrt(x * x + y * y);
}

void redshift(double r, double th, double ktkp, double& gg)
{
  double Omega;
    double uet;
    double g00,g03,g33,rderg00,rderg03,rderg33;
    double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42, t43, t44, t45, t46, t47, t48;
    double spin = chi;
    t1 = cos(th);
    t2 = r * r;
    t3 = pow(t2, 0.2e1);
    t4 = t2 * t3;
    t5 = r * t3;
    t6 = r * t2;
    t7 = spin * spin;
    t8 = t7 * pow(t1, 0.2e1);
    t9 = (t8 + t2) * r + eps3;
    t10 = a22 + t2;
    t11 = sin(th);
    t12 = t7 + t2;
    t13 = pow(t10, 0.2e1);
    t14 = pow(t11, 0.2e1);
    t15 = t7 * t13 * t14;
    t17 = a13 + t6;
    t18 = t7 * r;
    t19 = t18 * t10 * t14;
    t20 = t12 * t17;
    t21 = t20 - t19;
    t22 = 2.0 * r;
    t23 = -t22 + t12;
    t17 = pow(t12, 0.2e1) * pow(t17, 0.2e1) - t7 * t4 * t23 * t14;
    t24 = a22 + t7;
    t25 = r * a22;
    t26 = t7 * a13;
    t4 = 2.0 * t4 + t2 * (a13 * t24 + (a22 * t7 + (a13 + t25) * r) * r) + t26 * a22;
    t27 = a52 + t2;
    t28 = t23 * t3 - t15;
    t29 = t8 / 0.3e1 + t2;
    t30 = 0.2e1 / 0.3e1 * t7;
    t31 = 0.3e1 / 0.5e1 * t7;
    t32 = (t31 + t2) * r + 0.2e1 / 0.5e1 * a13;
    t31 = r * t32 - t31 * (a22 / 0.3e1 + t2) * t14;
    t33 = 0.1e1 / t21;
    t35 = pow(t33, 0.2e1);
    t36 = t33 * t35;
    t37 = 0.3e1 * t29;
    t38 = t35 * t14;
    t39 = t7 * a52;
    t40 = a52 + t7;
    t42 = 0.1e1 / r;
    t43 = t9 * t42;
    t44 = -2.0 * t11 * t35 * t1 * (-t43 * t17 + (t5 * t9 * t23 + t17) * t14 * t7) + 4.0 * t1 * t7 * t10 * t11 * t14 * t9 * t17 * t36;
    t21 = 0.1e1 / t21;
    t45 = 0.1e1 / t9;
    t46 = 0.1e1 / t27;
    t47 = 0.1e1 / t23;
    t48 = pow(t21, 0.2e1);
    t19 = 2.0 * t4 * t1 * spin * t11 * (t18 * t14 * (t19 + a22 * eps3 - ((-(a22 - t7) * r + a13 - eps3) * r - t8 * t10) * r - t26) + t20 * t9) * t21 * t48;
    t21 = 2.0 * t7 * t1;
    t5 = r * t9 * (-t12 * t3 + 2.0 * t5 + t15) * t48;
    t6 = (-2.0 * (t3 * (-t40 + r) + t8 * ((t2 * (r - 0.1e1) + a52) * r - t39)) * r + t3 * (a52 * (-6.0) - 0.3e1 * eps3) + (-t2 * t40 + t39) * eps3 + 4.0 * (eps3 + t39) * t6) * pow(t47, 0.2e1) * pow(t46, 0.2e1);

    g00 = t5;
    g03 = -t4 * t9 * spin * t14 * t48;
    g33 = t43 * t17 * t14 * t48;
    rderg00 = t35 * (t28 * (t9 * (0.10e2 * r * t31 * t33 - 0.1e1) - t37 * r) + (-6.0) * (t2 * ((-0.5e1 / 0.3e1 + r) * r + t30) - t30 * t10 * t14) * t2 * t9);
    rderg03 = t38 * spin * ((0.20e2 * t9 * t31 * t33 + t29 * (-6.0)) * t4 / 0.2e1 - 0.12e2 * r * ((t24 / 0.6e1 + t2 / 0.3e1) * a13 + t3 + t25 * (0.5e1 / 0.12e2 * t2 + t7 / 0.4e1)) * t9);
    rderg33 = 0.10e2 * t38 * t9 * (-t31 * t17 * t33 * t42 + t20 * t32 - 0.4e1 / 0.5e1 * t3 * ((-0.7e1 / 0.4e1 + r) * r + 0.3e1 / 0.4e1 * t7) * t7 * t14) + t38 * t17 * t42 * (-t43 + t37);
    
    Omega  = (-rderg03 + sqrt(rderg03*rderg03 - rderg00*rderg33))/rderg33;
    
    uet = sqrt(-g00 - 2.*g03*Omega - g33*Omega*Omega);
    
    gg = uet/(1. - ktkp*Omega);
}

double specific_energy(double r)
{
    double Omega, SE;
    double g00,g03,g33,rderg00,rderg03,rderg33;
    double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42, t43, t44, t45, t46, t47;
    double spin = chi;

    t1 = r * r;
    t2 = pow(t1, 0.2e1);
    t3 = pow(t2, 0.2e1);
    t4 = r * t1;
    t5 = r * t3;
    t6 = t1 * t2;
    t7 = r * t2;
    t8 = eps3 + t4;
    t9 = spin * spin;
    t10 = pow(t9, 0.2e1);
    t11 = t9 * t10;
    t12 = t9 * a22;
    t13 = a22 * a22;
    t14 = t9 * t13;
    t15 = t9 * a13;
    t16 = (r * (a13 + t4) - t12) * r + t15;
    t17 = r * a22;
    t18 = t15 * a22;
    t19 = a13 * a13;
    t20 = t10 * t19;
    t21 = t10 * a13;
    t22 = t19 * t9;
    t23 = eps3 / 0.2e1;
    t24 = a13 - t23;
    t25 = -0.10e2 * a13 + 0.8e1 * eps3;
    t26 = -0.18e2 * a13;
    t27 = 0.3e1 * eps3;
    t28 = a13 * eps3;
    t29 = t28 * t13;
    t30 = r * eps3;
    t31 = -0.12e2;
    t32 = 0.2e1 * eps3;
    t33 = 0.5e1 / 0.4e1 * t9;
    t34 = a22 / 0.4e1;
    t35 = 0.2e1 * t9;
    t36 = -a22 + t9;
    t37 = 0.4e1 * eps3;
    t38 = -t36;
    t39 = a13 + eps3;
    t40 = t22 * t30;
    t41 = 0.6e1 * a13;
    t42 = t9 + a13;
    t43 = t10 * a22;
    t44 = a13 + t32;
    t45 = eps3 / 0.6e1;
    t16 = 0.1e1 / t16;
    t46 = 0.1e1 / r;
    t47 = pow(t16, 0.2e1);
    t16 = t16 * t47;
    t33 = 0.4e1 * spin * (0.3e1 / 0.4e1 * t1 * (t4 * (a13 * t6 - t9 * (a13 * (a22 * t36 + t37) - t13 * eps3)) - t21 * a22 * t39) + t2 * t3 * (a22 + r) / 0.2e1 - t29 * t10 / 0.2e1 + t40 * t38 / 0.2e1 + (r * (t1 * (r * (r * (r * (r * (-0.5e1 / 0.2e1 * a13 + t12 + t32) + a13 * (t34 + t33) + 0.7e1 / 0.2e1 * a22 * (t9 + 0.5e1 / 0.14e2 * eps3)) + a13 * (-0.3e1 / 0.4e1 * a13 + 0.3e1 / 0.2e1 * eps3 - 0.9e1 / 0.2e1 * t9) + 0.3e1 / 0.2e1 * t14) + a13 * (-eps3 + t12 / 0.2e1) + 0.7e1 / 0.4e1 * t12 * eps3) + 0.2e1 * t28 * (0.7e1 / 0.8e1 * a22 + t9) + t19 * (-t34 - t35) + t12 * (t32 + t12)) + a13 * (-a13 * t9 * (a22 + t33) + t32 * t12)) + t13 * t10 * t39 / 0.4e1 + t23 * t38 * t19) * t4) * t16;
    t11 = (0.3e1 * t40 * (a13 * (r * t9 + t4) - t43) + t5 * (t4 * (t4 * (eps3 - t41 + t35) + a13 * (-t41 + t27) + t9 * (0.8e1 * t12 + t25)) + a13 * (a13 * (-0.2e1 * a13 + t27) - t37 * t9) + 0.20e2 * t43 * (a13 + eps3 / 0.4e1)) + t11 * a13 * t19 * eps3 - 0.10e2 * t6 * (a13 * (a13 * (-t28 / 0.10e2 + t11) - 0.2e1 / 0.5e1 * t43 * t44) - t12 * t3) - 0.6e1 * t9 * t7 * (t1 * (-a13 * (-t32 * t9 + a13 * (-a13 + eps3) + t43) + t24 * t6) + t19 * (t9 * (a13 + t45) - t45 * a22)) - 0.2e1 * t4 * (t2 * t3 * t4 + t20 * (a22 * eps3 + t44 * t9)) + t9 * t3 * (t1 * (t26 * (t42 - 0.2e1 / 0.3e1 * eps3) + 0.14e2 * t17 * (t23 + t42)) - 0.22e2 * t9 * (eps3 * (-0.4e1 / 0.11e2 * a22 - 0.9e1 / 0.22e2 * a13) + t19) + 0.4e1 * a13 * a22 * t44)) * t16 * pow(t46, 0.2e1);

    g00 = r * t8 * (0.2e1 * t1 * (t12 + t4) + t14 - t6) * t47;
    g03 = -spin * t8 * (0.2e1 * t6 + t1 * (a13 * (a22 + t9) + ((a13 + t17) * r + t12) * r) + t18) * t47;
    g33 = t8 * (t2 * (t2 * (t1 + t9) + t19) + t20 + 0.2e1 * t1 * ((t2 * (t9 + a13) + t21) * r + t22) + 0.4e1 * t15 * t7) * t47 * t46;
    rderg00 = -(0.2e1 * t2 * (t5 + t12 * (t14 + t28)) + t10 * (-t30 * a22 * t13 - t29) + t3 * (t1 * (0.6e1 * r * t24 + t25) + t9 * (0.6e1 * t13 + t26) + t27 * a13) - 0.4e1 * t4 * (t2 * (t12 * (a13 - 0.9e1 / 0.4e1 * eps3) + t28) + t24 * t13 * t10) + t9 * t1 * (t4 * (r * (0.8e1 * (0.7e1 / 0.8e1 * a13 + a22) * eps3 + 0.8e1 * t14 + (0.10e2 * a13 + 0.14e2 * a22) * t4) + (eps3 * (a13 - 0.3e1 / 0.4e1 * t13) + t18) * t31) - 0.6e1 * t28 * a22 * (-a22 / 0.2e1 + t9))) * t16;
    rderg03 = t33;
    rderg33 = -t11;
    
    Omega  = (-rderg03 + sqrt(rderg03*rderg03 - rderg00*rderg33))/rderg33;
    
    SE = -1.0*(g00 + Omega*g03)/sqrt(-1.0*g00 - 2.0*Omega*g03 - Omega*Omega*g33);
    
    return SE;
}

