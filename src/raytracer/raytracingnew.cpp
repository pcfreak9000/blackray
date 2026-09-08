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

static constexpr Real dobs = 1.0e+8; /* distance of the observer */

inline bool triviallyEnd(const PhaseVec &pvec, const int &iter,
    int &stop_integration) {
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

void raytrace(const size_t& ray_index, Record &rec,
    Env *env) {
  std::unique_ptr<InitialCondition> ic = initialcondition(ray_index);
  PhaseVec pvec = ic->pvec0;
  const Real b = ic->b;
  int &stop_integration = rec.stop_integration_condition;
  stop_integration = ST_INT_CONTINUE;
  PhaseVec pvecau;
  Real h = -1.0;
  Real prevh = -1.0;
  unsigned int iter = 0;
  bool freeze_h = false;

  /* ----- solve geodesic equations ----- */
  do {
    iter++;

    pvecau = pvec;
    adaptiveRK45<5>(b, pvec.vars, h, freeze_h, diffeqs);

    if (triviallyEnd(pvec, iter, stop_integration)) break;

    Entity *hit_ent = nullptr;
    bool hitb = env->checkIntersect(pvec.r, pvec.th, pvecau.r, pvecau.th,
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
      int nres = hit_ent->checkIntersect(pvec.r, pvec.th, pvecau.r, pvecau.th);
      if (!nres) {
        freeze_h = false;
        h = prevh;
        continue;
      }
      RayHit hit;
      hit.entity = hit_ent;
      hit.pvec = pvec;
      hit.pvecau = pvecau;
      //stop integration and "hit index"...
      stop_integration = nres;
      getRecord(rec, hit, ic.get());
    }
  } while (stop_integration == ST_INT_CONTINUE);
  rec.r = pvec.r;
  postRecord(rec, ic.get());
}
