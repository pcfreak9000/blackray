#include "raytracingnew.hpp"

#include <iostream>
#include <cmath>

#include "quadtree.hpp"
#include "def.hpp"
#include "diffeqs.hpp"
#include "redshift.hpp"
#include "metric.hpp"

#include "environment.hpp"

//in environment.cpp, clean up
void scalarProduct(Real met[4][4], Real *fvec0, Real *fvec1, Real &scal);
//void correct4VelNorm(Real met[4][4], Real norm, Real *fvel);

void raytrace(Real xobs, Real yobs, Real iobs, Real rin,
    Real disk_length_combined, RayHit &hit, int &stop_integration,
    Env* env) {
  Real dobs;
  Real xobs2, yobs2;
  Real hstart;
  Real r0, th0, phi0;
  Real kt0, kr0, kth0, kphi0;
  Real r02, s0, s02;
  Real fact1, fact2, fact3;
  Real r, th, phi;
  Real kr, kth;
  Real rau, thau, phiau, krau, kthau;
  Real const1;
  Real h;
  Real Delta;
  Real spin2 = spin * spin;

  Real carter, c02;
  Real b;

  Real met[4][4];
  Real diffs[5], vars[5], vars_temp[5], vars_4th[5], vars_5th[5], k1[5], k2[5],
      k3[5], k4[5], k5[5], k6[5];
  //Real xem[4];
  //Real gfactor;
  Real err, errmin, errmax;

  int check, check2 = 0;
  int i;

  /* ----- Set computational parameters ----- */
  dobs = 1.0e+8; /* distance of the observer */
  errmin = 1.0e-8;
  errmax = 1.0e-6;
  //atol = 1.0e-10;
  //rtol = 1.0e-10;
  //rtol = 1.0e3;
  Real thtol = 1.0e-8;
  int count, iter;

  hstart = -1.0;

  /* ----- compute photon initial conditions ----- */
  xobs2 = xobs * xobs;
  yobs2 = yobs * yobs;

  fact1 = yobs * std::sin(iobs) + dobs * std::cos(iobs);
  fact2 = dobs * std::sin(iobs) - yobs * std::cos(iobs);

  r02 = xobs2 + yobs2 + dobs * dobs;

  r0 = std::sqrt(r02);
  th0 = std::acos(fact1 / r0);
  phi0 = std::atan2(xobs, fact2);

  s0 = std::sin(th0);
  s02 = s0 * s0;

  kr0 = dobs / r0;
  kth0 = -(std::cos(iobs) - dobs * fact1 / r02)
      / std::sqrt(r02 - fact1 * fact1);
  kphi0 = -xobs * std::sin(iobs) / (xobs2 + fact2 * fact2);

  metric(r0, th0, met);

  fact3 = std::sqrt(
      met[0][3] * met[0][3] * kphi0 * kphi0
          - met[0][0]
              * (met[1][1] * kr0 * kr0 + met[2][2] * kth0 * kth0
                  + met[3][3] * kphi0 * kphi0));

  kt0 = -(met[0][3] * kphi0 + fact3) / met[0][0];

  b = -(met[3][3] * kphi0 + met[0][3] * kt0)
      / (met[0][0] * kt0 + met[0][3] * kphi0);

  kr0 /= fact3;
  kth0 /= fact3;

  /* ----- carter constant ----- */

  c02 = 1. - s02;

  carter = yobs2 - spin2 * c02 + xobs2 * c02;
  carter = std::sqrt(carter);

  /* ----- solve geodesic equations ----- */

  r = r0;
  th = th0;
  phi = phi0;

  kr = kr0;
  kth = kth0;

  //const0 = kt0;
  const1 = r02 * s02 * kphi0 / kt0;

  Real obsuarray[4] = { 1.0, 0.0, 0.0, 0.0 };
  Real obskarray[4] = { kt0, kr0, kth0, kphi0 };
  Real obsenergy;
  scalarProduct(met, obsuarray, obskarray, obsenergy);

  stop_integration = 0;

  h = hstart;
  count = 0;
  iter = 0;

  constexpr Real a1 = 1.0 / 4.0;
  constexpr Real b1 = 3.0 / 32.0;
  constexpr Real b2 = 9.0 / 32.0;
  constexpr Real c1 = 1932.0 / 2197.0;
  constexpr Real c2 = -7200.0 / 2197.0;
  constexpr Real c3 = 7296.0 / 2197.0;
  constexpr Real d1 = 439.0 / 216.0;
  constexpr Real d2 = -8.0;
  constexpr Real d3 = 3680.0 / 513.0;
  constexpr Real d4 = -845.0 / 4104.0;
  constexpr Real e1 = -8.0 / 27.0;
  constexpr Real e2 = 2.0;
  constexpr Real e3 = -3544.0 / 2565.0;
  constexpr Real e4 = 1859.0 / 4104.0;
  constexpr Real e5 = -11.0 / 40.0;
  constexpr Real f1 = 25.0 / 216.0;
  constexpr Real f2 = 0.0;
  constexpr Real f3 = 1408.0 / 2565.0;
  constexpr Real f4 = 2197.0 / 4104.0;
  constexpr Real f5 = -1.0 / 5.0;
  constexpr Real g1 = 16.0 / 135.0;
  constexpr Real g2 = 0.0;
  constexpr Real g3 = 6656.0 / 12825.0;
  constexpr Real g4 = 28561.0 / 56430.0;
  constexpr Real g5 = -9.0 / 50.0;
  constexpr Real g6 = 2.0 / 55.0;
  Real prevh = -1.0;

  do {
    iter++;
    vars[0] = r;
    vars[1] = th;
    vars[2] = phi;
    vars[3] = kr;
    vars[4] = kth;
//evolve the system one step such that the error stays below certain limits. For this, adaptively increase or decrease step size
    do {
      check = 0;

      /* ----- compute RK1 ----- */

      diffeqs(b, vars, diffs);
      for (i = 0; i <= 4; i++) {
        k1[i] = h * diffs[i];
        vars_temp[i] = vars[i] + a1 * k1[i];
      }

      /* ----- compute RK2 ----- */

      diffeqs(b, vars_temp, diffs);
      for (i = 0; i <= 4; i++) {
        k2[i] = h * diffs[i];
        vars_temp[i] = vars[i] + b1 * k1[i] + b2 * k2[i];
      }

      /* ----- compute RK3 ----- */

      diffeqs(b, vars_temp, diffs);
      for (i = 0; i <= 4; i++) {
        k3[i] = h * diffs[i];
        vars_temp[i] = vars[i] + c1 * k1[i] + c2 * k2[i] + c3 * k3[i];
      }

      /* ----- compute RK4 ----- */

      diffeqs(b, vars_temp, diffs);
      for (i = 0; i <= 4; i++) {
        k4[i] = h * diffs[i];
        vars_temp[i] = vars[i] + d1 * k1[i] + d2 * k2[i] + d3 * k3[i]
            + d4 * k4[i];
      }

      /* ----- compute RK5 ----- */

      diffeqs(b, vars_temp, diffs);
      for (i = 0; i <= 4; i++) {
        k5[i] = h * diffs[i];
        vars_temp[i] = vars[i] + e1 * k1[i] + e2 * k2[i] + e3 * k3[i]
            + e4 * k4[i] + e5 * k5[i];
      }

      /* ----- compute RK6 ----- */

      diffeqs(b, vars_temp, diffs);
      for (i = 0; i <= 4; i++)
        k6[i] = h * diffs[i];

      /* ----- local error ----- */

      for (i = 0; i <= 4; i++) {
        vars_4th[i] = vars[i] + f1 * k1[i] + f2 * k2[i] + f3 * k3[i]
            + f4 * k4[i] + f5 * k5[i];
        vars_5th[i] = vars[i] + g1 * k1[i] + g2 * k2[i] + g3 * k3[i]
            + g4 * k4[i] + g5 * k5[i] + g6 * k6[i];

        err = fabs(
            (vars_4th[i] - vars_5th[i]) / std::max(vars_4th[i], vars[i]));

        if (err > errmax && check2 == 0) check = 1;
        else if (err < errmin && check != 1 && check2 == 0) check = -1;
      }

      if (check == 1) {
        h /= 2.0;
#ifdef ITER_WARN
        if (iter > MAX_ITER - 10) {
          std::cout << "descale0" << std::endl;
        }
#endif
      } else if (check == -1) h *= 2.0;

    } while (check == 1);

    /* ----- solutions to the fourth-order RKN method ----- */
//apply the new step to the variables
    rau = r;
    thau = th;
    phiau = phi;
    krau = kr;
    kthau = kth;

    r = vars_4th[0];
    th = vars_4th[1];
    phi = vars_4th[2];
    kr = vars_4th[3];
    kth = vars_4th[4];
    Delta = SQR(r) - 2.0 * r + spin2;
//check if the new position ends the integration
    if (Delta < 1.0e-3) {
      stop_integration = 4; // printf("photon crosses the horizon\n"); /* the
      // photon crosses the horizon */
      break;
    }
    if (r < 1.0) {
      stop_integration = 5; // printf("photon crosses the horizon\n"); /* the
      // photon crosses the horizon */
      break;
    }

    if (r != r) {
      stop_integration = 6; // printf("numerical problem\n");          /*
      // numerical problems! */
      break;
    }

    if (r > 1.05 * dobs) {
      stop_integration = 7; // printf("photon escaped to infinity\n");   /* the
      // photon escapes to infinity */
      break;
    }
#ifdef ITER_WARN
    if (iter > MAX_ITER - 10) {
      std::cout << "Reaching max iter with..." << std::endl;
      std::cout << h << " " << iter << std::endl;
      std::cout << r << " " << th << " " << phi << std::endl;
      std::cout << rau << " " << thau << " " << phiau << std::endl;
    }
#endif
    if (iter > MAX_ITER) {
      stop_integration = 255;
      break;
    }
    int res = NO_INTERSECT;
    SurfacePoint spi;
    Entity *hit_ent = env->checkIntersect(r, th, rau, thau, spi, res);

    //deal with (possible) intersection
    if (res == INTERSECT) {
#ifdef ITER_WARN
      if (iter > MAX_ITER - 10) {
        std::cout << "int" << std::endl;
      }
#endif
      if (check2 != 1) {
        prevh = h;
      }
      check2 = 1; //don't adapt stepsize anymore, this is now done manually to reach certain tolerances
      if (std::fabs(th - thau) <= thtol) {
        count++;
      }
      if (count > 0) {
        int nres = hit_ent->checkIntersect(r, th, rau, thau, spi);
        if (nres==INTERSECT) {
          IntegratorData id; //oof
          id.b = b;
          id.carter = carter;
          id.const1 = const1;
          id.kr = kr;
          id.kth = kth;
          id.obsenergy = obsenergy;
          id.r = r;
          id.rprev = rau;
          id.th = th;
          id.thprev = thau;
          stop_integration = hit_ent->intersect(id, spi, hit);
          if (hit.gfactor < 0.0) {
            std::cout << "gfactor is < 0.0, ignoring ray" << std::endl;
            stop_integration = 6;
            hit.r = r;
            hit.gfactor = 1.0;
            hit.cosem = 0.0;
          }
          if(std::isnan(hit.gfactor)) {
            std::cout << "gfactor is nan, ignoring ray" << std::endl;
            stop_integration = 6;
            hit.r = r;
            hit.gfactor = 1.0;
            hit.cosem = 0.0;
          }
#ifdef DEBUG_COSEM
          if (std::isnan(hit.cosem)) {
            std::cout << "Cosem is nan, ignoring ray" << std::endl;
            stop_integration = 6;
            hit.r = r;
            hit.gfactor = 1.0;
            hit.cosem = 0.0;
          }
          if (hit.cosem > 1.05) {
            std::cout << "Cosem > 1.05 detected, ignoring ray: " << hit.cosem
                << std::endl;
            stop_integration = 6;
            hit.r = r;
            hit.gfactor = 1.0;
            hit.cosem = 0.0;
          }
          if (hit.cosem > 1.0) {
            std::cout << "1.05 >= Cosem > 1.0 detected, clamping to 1.0: "
                << hit.cosem << std::endl;
            hit.cosem = 1.0;
          }
#endif
        } else {
          check2 = 0;
          count = 0;
          h = prevh;
          // this is a simplification, we can only be sure about the final
          // whereabouts of the photon if it ends up beyond the event horizon ore
          // is ejected to infinity
          //stop_integration = 2; /* the photon misses the disk */
        }
      } else {
        r = rau;
        th = thau;
        phi = phiau;
        kr = krau;
        kth = kthau;
        h /= 2.0;
#ifdef ITER_WARN
        if (iter > MAX_ITER - 10) {
          std::cout << "descale1" << std::endl;
        }
#endif
      }
    }
  } while (stop_integration == 0);
  if (stop_integration != 512
      && (stop_integration < 128 || stop_integration > 131) && stop_integration != 600) {
    hit.r = r;
    hit.gfactor = 1.0;
    hit.cosem = 0.0;
  }

}
