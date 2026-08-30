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

static constexpr Real a1 = 1.0 / 4.0;
static constexpr Real b1 = 3.0 / 32.0;
static constexpr Real b2 = 9.0 / 32.0;
static constexpr Real c1 = 1932.0 / 2197.0;
static constexpr Real c2 = -7200.0 / 2197.0;
static constexpr Real c3 = 7296.0 / 2197.0;
static constexpr Real d1 = 439.0 / 216.0;
static constexpr Real d2 = -8.0;
static constexpr Real d3 = 3680.0 / 513.0;
static constexpr Real d4 = -845.0 / 4104.0;
static constexpr Real e1 = -8.0 / 27.0;
static constexpr Real e2 = 2.0;
static constexpr Real e3 = -3544.0 / 2565.0;
static constexpr Real e4 = 1859.0 / 4104.0;
static constexpr Real e5 = -11.0 / 40.0;
static constexpr Real f1 = 25.0 / 216.0;
static constexpr Real f2 = 0.0;
static constexpr Real f3 = 1408.0 / 2565.0;
static constexpr Real f4 = 2197.0 / 4104.0;
static constexpr Real f5 = -1.0 / 5.0;
static constexpr Real g1 = 16.0 / 135.0;
static constexpr Real g2 = 0.0;
static constexpr Real g3 = 6656.0 / 12825.0;
static constexpr Real g4 = 28561.0 / 56430.0;
static constexpr Real g5 = -9.0 / 50.0;
static constexpr Real g6 = 2.0 / 55.0;

static constexpr Real errmin = 1.0e-8;
static constexpr Real errmax = 1.0e-6;

static constexpr Real dobs = 1.0e+8; /* distance of the observer */
//atol = 1.0e-10;
//rtol = 1.0e-10;
//rtol = 1.0e3;
static constexpr Real thtol = 1.0e-8;

union PhaseVec {
  Real vars[5];
  struct {
    Real r, th, phi;
    Real kr, kth;
  };
};

template<size_t N>
void adaptiveRK45(const Real &b, Real vars[N], Real &h, const bool &freeze_h) {
  //evolve the system one step such that the error stays below certain limits. For this, adaptively increase or decrease step size
  Real diffs[N], vars_4th[N], vars_temp[N], k1[N], k2[N], k3[N], k4[N], k5[N]; //, k6[5];
  do {
    int check = 0;

    /* ----- compute RK1 ----- */
    //maybe this diffeqs only needs computation once and not everytime to adjust h?
    diffeqs(b, vars, diffs);
#pragma omp simd
    for (size_t i = 0; i < N; i++) {
      k1[i] = h * diffs[i];
      vars_temp[i] = vars[i] + a1 * k1[i];
    }

    /* ----- compute RK2 ----- */

    diffeqs(b, vars_temp, diffs);
#pragma omp simd
    for (size_t i = 0; i < N; i++) {
      k2[i] = h * diffs[i];
      vars_temp[i] = vars[i] + b1 * k1[i] + b2 * k2[i];
    }

    /* ----- compute RK3 ----- */

    diffeqs(b, vars_temp, diffs);
#pragma omp simd
    for (size_t i = 0; i < N; i++) {
      k3[i] = h * diffs[i];
      vars_temp[i] = vars[i] + c1 * k1[i] + c2 * k2[i] + c3 * k3[i];
    }

    /* ----- compute RK4 ----- */

    diffeqs(b, vars_temp, diffs);
#pragma omp simd
    for (size_t i = 0; i < N; i++) {
      k4[i] = h * diffs[i];
      vars_temp[i] = vars[i] + d1 * k1[i] + d2 * k2[i] + d3 * k3[i]
          + d4 * k4[i];
    }

    /* ----- compute RK5 ----- */

    diffeqs(b, vars_temp, diffs);
#pragma omp simd
    for (size_t i = 0; i < N; i++) {
      k5[i] = h * diffs[i];
      vars_temp[i] = vars[i] + e1 * k1[i] + e2 * k2[i] + e3 * k3[i] + e4 * k4[i]
          + e5 * k5[i];
    }

    /* ----- compute RK6 ----- */

    diffeqs(b, vars_temp, diffs);
//this for loop is integrated into the next for-loop because it is redundant
//    for (int i = 0; i <= 4; i++)
//      k6[i] = h * diffs[i];

    /* ----- local error ----- */

    for (size_t i = 0; i < N; i++) {
      vars_4th[i] = vars[i] + f1 * k1[i] + f2 * k2[i] + f3 * k3[i] + f4 * k4[i]
          + f5 * k5[i];
      //vars_5th[i] =
      Real varfith = vars[i] + g1 * k1[i] + g2 * k2[i] + g3 * k3[i] + g4 * k4[i]
          + g5 * k5[i] + g6 * diffs[i] * h; //k6[i];

      Real err = std::fabs(
          (vars_4th[i] - varfith) / std::max(vars_4th[i], vars[i]));

      if (err > errmax && !freeze_h) check = 1;
      else if (err < errmin && check != 1 && !freeze_h) check = -1;
    }

    if (check == 1) {
      h /= 2.0;
    } else {
      if (check == -1) h *= 2.0;
      //apply the new step to the variables
      for (size_t i = 0; i < N; i++) {
        vars[i] = vars_4th[i];
      }
      break;
    }

  } while (true);
}
void invalidRay(const PhaseVec &pvec, RayHit& hit){
  hit.r = pvec.r;
  hit.gfactor = 1.0;
  hit.cosem = 0.0;
}
bool triviallyEnd(const PhaseVec &pvec, const int &iter,
    int &stop_integration) {
  //check if the new position ends the integration
  Real Delta = SQR(pvec.r) - 2.0 * pvec.r + SQR(spin);
  if (Delta < 1.0e-3) {
    stop_integration = 4; // printf("photon crosses the horizon\n"); /* the
    // photon crosses the horizon */
    return true;
  }
  if (pvec.r < 1.0) {
    stop_integration = 5; // printf("photon crosses the horizon\n"); /* the
    // photon crosses the horizon */
    return true;
  }

  if (std::isnan(pvec.r)) {
    stop_integration = 6; // printf("numerical problem\n");          /*
    // numerical problems! */
    return true;
  }

  if (pvec.r > 1.05 * dobs) {
    stop_integration = 7; // printf("photon escaped to infinity\n");   /* the
    // photon escapes to infinity */
    return true;
  }
  if (iter > MAX_ITER) {
    stop_integration = 255;
    return true;
  }
  return false;
}



void handleNonsense(const PhaseVec &pvec, RayHit &hit, int &stop_integration) {
  if (hit.gfactor < 0.0) {
    std::cout << "gfactor is < 0.0, ignoring ray" << std::endl;
    stop_integration = 6;
    invalidRay(pvec, hit);
  }
  if (std::isnan(hit.gfactor)) {
    std::cout << "gfactor is nan, ignoring ray" << std::endl;
    stop_integration = 6;
    invalidRay(pvec, hit);
  }
#ifdef DEBUG_COSEM
  if (std::isnan(hit.cosem)) {
    std::cout << "Cosem is nan, ignoring ray" << std::endl;
    stop_integration = 6;
    invalidRay(pvec, hit);
  }
  if (hit.cosem > 1.05) {
    std::cout << "Cosem > 1.05 detected, ignoring ray: " << hit.cosem
        << std::endl;
    stop_integration = 6;
    invalidRay(pvec, hit);
  }
  if (hit.cosem > 1.0) {
    std::cout << "1.05 >= Cosem > 1.0 detected, clamping to 1.0: " << hit.cosem
        << std::endl;
    hit.cosem = 1.0;
  }
#endif
}

void raytrace(Real xobs, Real yobs, Real iobs, Real rin,
    Real disk_length_combined, RayHit &hit, int &stop_integration, Env *env) {
  /* ----- compute photon initial conditions ----- */
  const Real xobs2 = xobs * xobs;
  const Real yobs2 = yobs * yobs;

  const Real fact1 = yobs * std::sin(iobs) + dobs * std::cos(iobs);
  const Real fact2 = dobs * std::sin(iobs) - yobs * std::cos(iobs);

  const Real r02 = xobs2 + yobs2 + dobs * dobs;

  const Real r0 = std::sqrt(r02);
  const Real th0 = std::acos(fact1 / r0);
  const Real phi0 = std::atan2(xobs, fact2);

  const Real s0 = std::sin(th0);
  const Real s02 = s0 * s0;

  const Real kr0_unscl = dobs / r0;
  const Real kth0_unscl = -(std::cos(iobs) - dobs * fact1 / r02)
      / std::sqrt(r02 - fact1 * fact1);
  const Real kphi0 = -xobs * std::sin(iobs) / (xobs2 + fact2 * fact2);

  Real met[4][4];//kind of a waste of space
  metric(r0, th0, met);

  const Real fact3 = std::sqrt(
      met[0][3] * met[0][3] * kphi0 * kphi0
          - met[0][0]
              * (met[1][1] * kr0_unscl * kr0_unscl + met[2][2] * kth0_unscl * kth0_unscl
                  + met[3][3] * kphi0 * kphi0));

  const Real kt0 = -(met[0][3] * kphi0 + fact3) / met[0][0];

  const Real b = -(met[3][3] * kphi0 + met[0][3] * kt0)
      / (met[0][0] * kt0 + met[0][3] * kphi0);

  const Real kr0 = kr0_unscl / fact3;
  const Real kth0 = kth0_unscl / fact3;

  /* ----- some more constants ----- */

  const Real c02 = 1. - s02;
  const Real carter = std::sqrt(yobs2 - SQR(spin) * c02 + xobs2 * c02);
  //const0 = kt0;
  const Real const1 = r02 * s02 * kphi0 / kt0;

  /* ----- prepare solver ----- */

  PhaseVec pvec, pvecau;
  pvec.r = r0;
  pvec.th = th0;
  pvec.phi = phi0;

  pvec.kr = kr0;
  pvec.kth = kth0;

  Real obsuarray[4] = { 1.0, 0.0, 0.0, 0.0 };
  Real obskarray[4] = { kt0, kr0, kth0, kphi0 };
  Real obsenergy;
  scalarProduct(met, obsuarray, obskarray, obsenergy);

  stop_integration = 0;

  Real h = -1.0;
  Real prevh = -1.0;
  unsigned int iter = 0;
  bool freeze_h = false;

  /* ----- solve geodesic equations ----- */
  do {
    iter++;

    pvecau = pvec;
    adaptiveRK45<5>(b, pvec.vars, h, freeze_h);

    if (triviallyEnd(pvec, iter, stop_integration)){
      invalidRay(pvec, hit);
      break;
    }

    SurfacePoint spi;
    Entity *hit_ent = nullptr;
    bool hitb = env->checkIntersect(pvec.r, pvec.th, pvecau.r, pvecau.th, spi,
        hit_ent);

    if (hitb) {
      //don't adapt stepsize anymore, this is now done manually to reach certain tolerances
      if (!freeze_h) {
        prevh = h;
        freeze_h = true;
      }
      if (std::fabs(pvec.th - pvecau.th) > thtol) {
        pvec = pvecau;
        h /= 2.0;
        continue;
      }
      bool nres = hit_ent->checkIntersect(pvec.r, pvec.th, pvecau.r, pvecau.th,
          spi);
      if (!nres) {
        freeze_h = false;
        h = prevh;
        continue;
      }
      IntegratorData id; //oof
      id.b = b;
      id.carter = carter;
      id.const1 = const1;
      id.kr = pvec.kr;
      id.kth = pvec.kth;
      id.obsenergy = obsenergy;
      id.r = pvec.r;
      id.rprev = pvecau.r;
      id.th = pvec.th;
      id.thprev = pvecau.th;
      stop_integration = hit_ent->intersect(id, spi, hit);
      handleNonsense(pvec, hit, stop_integration);
    }
  } while (stop_integration == 0);

}
