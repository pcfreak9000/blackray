#include <iostream>
#include <cmath>

#include "raytracingnew.hpp"

#include "aRK45.hpp"
#include "diffeqs.hpp"

#include "problem.hpp"

//atol = 1.0e-10;
//rtol = 1.0e-10;
//rtol = 1.0e3;
static constexpr Real thtol = 1.0e-8;

inline bool triviallyEnd(const PhaseVec &pvec, const int &iter,
    int &stop_integration, const Real& dobs) {
  //check if the new position ends the integration
  Real Delta = SQR(pvec.r) - 2.0 * pvec.r + SQR(spin);
  if (Delta < 1.0e-3) {
    stop_integration = ST_INT_PHOTON_HORIZON_CROSS_DCRIT; // printf("photon crosses the horizon\n"); /* the
    // photon crosses the horizon */
    return true;
  }
  if (pvec.r < 1.0) {
    stop_integration = ST_INT_PHOTON_HORIZON_CROSS_RCRIT; // printf("photon crosses the horizon\n"); /* the
    // photon crosses the horizon */
    return true;
  }

  if (std::isnan(pvec.r)) {
    stop_integration = ST_INT_PROBLEM; // printf("numerical problem\n");          /*
    // numerical problems! */
    return true;
  }

  if (pvec.r > 1.05 * dobs) {
    stop_integration = ST_INT_PHOTON_ESCAPES; // printf("photon escaped to infinity\n");   /* the
    // photon escapes to infinity */
    return true;
  }
  if (iter > MAX_ITER) {
    stop_integration = ST_INT_MAX_ITERATIONS;
    return true;
  }
  return false;
}

void raytrace(const size_t &ray_index, Env *env) {
  std::unique_ptr<InitialCondition> ic = initialcondition(ray_index);
  PhaseVec pvec = ic->pvec0;
  const Real b = ic->b;
  int stop_integration = ST_INT_CONTINUE;
  PhaseVec pvecau;
  Real h = -1.0;
  Real prevh = -1.0;
  unsigned int iter = 0;
  bool freeze_h = false;
  RayHit hit;
  hit.entity = nullptr;
  /* ----- solve geodesic equations ----- */
  do {
    iter++;

    pvecau = pvec;
    adaptiveRK45<5>(b, pvec.vars, h, freeze_h, diffeqs);

    if (triviallyEnd(pvec, iter, stop_integration, ic->dobs)) break;

    Entity *hit_ent = nullptr;
    bool hitb = env->checkIntersect(pvec.r, pvec.th, pvecau.r, pvecau.th,
        hit_ent);

    if (!hitb) continue;
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
    int nres = hit_ent->checkIntersect(pvec.r, pvec.th, pvecau.r, pvecau.th);
    if (!nres) {
      freeze_h = false;
      h = prevh;
      continue;
    }
    hit.entity = hit_ent;
    stop_integration = nres;
  } while (stop_integration == ST_INT_CONTINUE);
  hit.pvec = pvec;
  hit.pvecau = pvecau;
  notifyDone(stop_integration, ic.get(), hit, ray_index);
}

void raytrace(Env *env, size_t ray_count_total) {
  std::cout << "Starting raytracing loop" << std::endl;
#pragma omp parallel for ordered schedule(dynamic)
  for (size_t ray_index = 0; ray_index < ray_count_total; ray_index++) {
    /* ----- unfold initial conditions ----- */
    if (ray_index % 10000 == 0) {
      std::cout << "Progress: " << ray_index / (double) (ray_count_total)
          << std::endl;
    }
    /* ----- raytrace ----- */
    raytrace(ray_index, env);
  }
}

