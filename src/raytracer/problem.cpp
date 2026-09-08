#include <cmath>
#include <iostream>
#include <memory>

#include "problem.hpp"

#include "raytracingnew.hpp"
#include "metric.hpp"
#include "find_isco.hpp"
#include "quadtree.hpp"
#include "environment.hpp"

//in environment.cpp, clean up
void scalarProduct(Real met[4][4], Real *fvec0, Real *fvec1, Real &scal);
static constexpr Real dobs = 1.0e+8; /* distance of the observer double def oof */
static myinitialconditions initcons;
void initialConditionGenerator(const Real &rstep, const Real &pstep,
    myinitialconditions &ic) {
  /* ----- Set computational parameters ----- */

  const Real robs_i = 0.1; //this was 1, but that is too big and will yield artifacts
  const Real robs_f = 215;

  // rstep  = 1.008;
  //const Real rstep2 = (rstep - 1) / rstep;
  // pstep  = 2*Pi/720;

  for (Real pobs = 0; pobs < 2 * Pi - 0.5 * pstep; pobs = pobs + pstep) {
    ic.pobs_vec.push_back(pobs);
  }

  for (Real robs = robs_i; robs < robs_f; robs = robs * rstep) {
    ic.robs_vec.push_back(robs);
  }
  ic.count_total = ic.robs_vec.size() * ic.pobs_vec.size();

}
void setupProblem(int argc, char *argv[], Env *env, size_t &ray_count_total) {
  /* ----- SETUP ENVIRONMENT ----- */
  Real maxx;
  Real maxy;
  const char *diskdatafile = argv[10];
  const Real rstep = atof(argv[8]);
  const Real pstep = atof(argv[9]);
  const Real spin2 = SQR(spin);
  std::unique_ptr<QuadTree> tree = readFileToTree(diskdatafile, maxx, maxy);
//  if (!tree) return 1;
  Real maxr_xdir = std::sqrt(SQR(maxx+10.0) - spin2) + 10.0;
  Real maxr_ydir = std::sqrt(SQR(maxy + 10.0)) + 10.0;
  Real checkr = maxr_ydir;
  if (maxr_xdir > maxr_ydir) checkr = maxr_xdir;
  Real isco;

  find_isco(15.0, isco); /* Depends upon the properties of BH */
//  GRMHDDisk disk(tree.get(), checkr);
  //env.addEntity(&disk);
  env->addEntity(std::make_unique<ThinDisk>(isco, 200));
  env->addEntity(std::make_unique<PlungingRegion>(isco));

  initialConditionGenerator(rstep, pstep, initcons);
  ray_count_total = initcons.count_total;
}

std::unique_ptr<InitialCondition> initialcondition(const size_t &ray_index) {
  std::unique_ptr<gf_ic> cond = std::make_unique<gf_ic>();
  /* ----- compute photon initial conditions ----- */
  const Real iobs = Pi / 180 * iobs_deg; /* inclination angle of the observer in rad */
  const size_t robs_index = ray_index / initcons.pobs_vec.size();
  const size_t pobs_index = ray_index % initcons.pobs_vec.size();
  const Real robs = initcons.robs_vec[robs_index];
  const Real pobs = initcons.pobs_vec[pobs_index];
  cond->robs = robs;
  cond->pobs = pobs;

  const Real xobs = robs * std::cos(pobs);
  const Real yobs = robs * std::sin(pobs);

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

  Real met[4][4];    //kind of a waste of space
  metric(r0, th0, met);

  const Real fact3 = std::sqrt(
      met[0][3] * met[0][3] * kphi0 * kphi0
          - met[0][0]
              * (met[1][1] * kr0_unscl * kr0_unscl
                  + met[2][2] * kth0_unscl * kth0_unscl
                  + met[3][3] * kphi0 * kphi0));

  const Real kt0 = -(met[0][3] * kphi0 + fact3) / met[0][0];

  cond->b = -(met[3][3] * kphi0 + met[0][3] * kt0)
      / (met[0][0] * kt0 + met[0][3] * kphi0);

  const Real kr0 = kr0_unscl / fact3;
  const Real kth0 = kth0_unscl / fact3;

  /* ----- some more constants ----- */

  const Real c02 = 1. - s02;
  const Real carter = std::sqrt(yobs2 - SQR(spin) * c02 + xobs2 * c02);
  //const0 = kt0;
  const Real const1 = r02 * s02 * kphi0 / kt0;
  cond->carter = carter;
  cond->const1 = const1;
  /* ----- prepare solver ----- */

  cond->pvec0.r = r0;
  cond->pvec0.th = th0;
  cond->pvec0.phi = phi0;

  cond->pvec0.kr = kr0;
  cond->pvec0.kth = kth0;

  Real obsuarray[4] = { 1.0, 0.0, 0.0, 0.0 };
  Real obskarray[4] = { kt0, kr0, kth0, kphi0 };
  scalarProduct(met, obsuarray, obskarray, cond->obsenergy);

  return cond;
}
inline void invalidRay(Record &rec) {
  rec.cosem = 0.0;
  rec.gfactor = 1.0;
  rec.stop_integration_condition = ST_INT_PROBLEM;
}

inline void handleNonsense(Record &rec) {
  if (rec.gfactor < 0.0) {
    std::cout << "gfactor is < 0.0, ignoring ray" << std::endl;
    invalidRay(rec);
  }
  if (std::isnan(rec.gfactor)) {
    std::cout << "gfactor is nan, ignoring ray" << std::endl;
    invalidRay(rec);
  }
#ifdef DEBUG_COSEM
  if (std::isnan(rec.cosem)) {
    std::cout << "Cosem is nan, ignoring ray" << std::endl;
    invalidRay(rec);
  }
  if (rec.cosem > 1.05) {
    std::cout << "Cosem > 1.05 detected, ignoring ray: " << rec.cosem
        << std::endl;
    invalidRay(rec);
  }
  if (rec.cosem > 1.0) {
    std::cout << "1.05 >= Cosem > 1.0 detected, clamping to 1.0: " << rec.cosem
        << std::endl;
    rec.cosem = 1.0;
  }
#endif
}

void postRecord(Record &rec, const InitialCondition *const ic) {
  const gf_ic *const gic = (gf_ic*) ic;
  rec.xobs = gic->robs * std::cos(gic->pobs);
  rec.yobs = gic->robs * std::sin(gic->pobs);
  rec.robs = gic->robs;
}

void getRecord(Record &rec, const RayHit &hit,
    const InitialCondition *const ic) {
  Real gfactor;
  Real cosem;
  int x = hit.entity->calculateRedshift(ic, hit, gfactor, cosem);
  if (x != 0) rec.stop_integration_condition = ST_INT_PROBLEM;
  rec.cosem = cosem;
  rec.gfactor = gfactor;
  handleNonsense(rec);
}
