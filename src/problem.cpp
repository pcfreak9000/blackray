#include <cmath>
#include <iostream>
#include <fstream>
#include <memory>

#include "problem.hpp"

#include "raytracingnew.hpp"
#include "metric.hpp"
#include "find_isco.hpp"
#include "quadtree.hpp"
#include "environment.hpp"

static constexpr Real dobs = 1.0e+8; /* distance of the observer */

struct Record {
  int stop_integration_condition;
  size_t photon_index;
  Real xobs, yobs;
  Real r;
  Real gfactor;
  Real cosem;
  bool output;
  size_t ray_index;
  Real robs;
};

static myinitialconditions initcons;
static std::vector<Record> records;
static size_t photon_index = 0;
static Real iobs_deg;


void ironlinestuff(const Real& rstep, const Real& alpha, const char *tempdirs, std::vector<Record> &records){

  /* ----- Set model for the spectral line ----- */

  constexpr Real E_line = 6.4; /* energy rest of the line in keV */
  constexpr Real N_0 = 1.0; /* normalization */
  // alpha  = -3;     radial power law index

  const Real rstep2 = (rstep - 1) / rstep;

  Real E_obs[IMAX];
  Real N_obs[IMAX];
  E_obs[0] = 0.0125000002; /* minimum photon energy detected by the observer; in keV */
  N_obs[0] = 0;
  constexpr Real E_step = 0.025;
  for (size_t i = 1; i < IMAX; i++) {
    E_obs[i] = E_obs[i - 1] + E_step;
    N_obs[i] = 0;
  }

  for (Record rec : records) {
    if (rec.output) {
      Real &gfactor = rec.gfactor;

      /* --- integration - part 1 --- */
      Real pp = gfactor * E_line;
      Real bucket = (pp - E_obs[0]) / E_step;
      if (bucket >= 0.0) {
        size_t index = std::floor(bucket);
        if (index < IMAX) {
          Real qq = gfactor * gfactor * gfactor * gfactor;
          qq = qq * std::pow(rec.r, alpha);
          /* --- integration - part 2 --- */
          //N_obs_add =
          N_obs[index] += SQR(rec.robs) * rstep2 * qq;
        }
      }
    }
  }

  Real N_tot = 0.0;
#pragma omp parallel for reduction(+:N_tot)
  for (size_t i = 0; i < IMAX; i++) {
    N_obs[i] = N_0 * N_obs[i] / E_obs[i];
    N_tot += N_obs[i];
  }

  /* --- print iron line --- */
  char filename_o[256];
  /*Iron line output file*/
  // sprintf(filename_o,"iron_a%.03f.epsilon_r%.02f.epsilon_t%.02f.i%.02f.dat",spin,epsi3,iobs_deg);
  // sprintf(filename_o,"ironline_data/iron_a%.05Le.i%.02Le.e_%.02Le.a13_%.02Le.a22_%.02Le.a52_%.02Le.dat",spin,iobs_deg,epsi3,a13,a22,a52);
  std::string s1("ironline_data/"
      "iron_a_%.05Lf_i_%.05Lf_e_%.05Lf_a13_%.05Lf_a22_%.05Lf_a52_%.05Lf.dat");
  snprintf(filename_o, sizeof(filename_o), (tempdirs + s1).c_str(), spin,
      iobs_deg, epsi3, a13, a22, a52);
  FILE *foutput = fopen(filename_o, "w");
  if (foutput == nullptr) std::cerr << "Problems with iron line file!"
      << std::endl;

  for (size_t i = 0; i < IMAX; i++) {
    fprintf(foutput, "%Lf %.10Lf\n", E_obs[i], N_obs[i] / N_tot);
  }
  fclose(foutput);
}

void write(const char *tempdir, const char *outtxt, std::vector<Record> &recs) {
  char filename_o2[256];

  FILE *foutput_coord;

  std::string tempdirs(tempdir);

  /*photon data output file*/
  // sprintf(filename_o2,"coord_a%.03f.epsilon_r%.02f.epsilon_t%.02f.i%.02f.dat",spin,epsi3,iobs_deg);
  std::string s2(
      "data/"
          "photons_data_a%.05Lf_i_%.05Lf_e_%.05Lf_a13_%.05Lf_a22_%.05Lf_a52_%.05Lf.dat");
  snprintf(filename_o2, sizeof(filename_o2), (tempdirs + s2).c_str(), spin,
      iobs_deg, epsi3, a13, a22, a52);

  foutput_coord = fopen(filename_o2, "w");
  if (foutput_coord == nullptr) std::cerr << "Problems with data file!"
      << std::endl;

  std::ofstream tmpOutFile(outtxt);
  for (Record r : recs) {
    if (r.output) {
      fprintf(foutput_coord, "%zu %Lf %Lf %Lf %Lf %Lf\n", r.photon_index,
          r.xobs, r.yobs, r.r, r.gfactor, r.cosem);
      if (!RESTRICT_DEBUGFILE_CRIT && r.ray_index % DEBUGFILE_OUT_DIV == 0) {
        tmpOutFile << r.xobs << " " << r.yobs << " " << r.gfactor << " "
            << r.stop_integration_condition << " " << std::endl;
      }
    } else {
      if ((!RESTRICT_DEBUGFILE_CRIT && r.ray_index % DEBUGFILE_OUT_DIV == 0)
          || (r.stop_integration_condition == 255
              || r.stop_integration_condition == 6)) {
        tmpOutFile << r.xobs << " " << r.yobs << " " << 1.0 << " "
            << r.stop_integration_condition << std::endl;
      }
    }
  }
  tmpOutFile.close();
  fclose(foutput_coord);
}

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
  Real maxy;static Real iobs_deg;

  const char *diskdatafile = argv[10];
  const Real rstep = atof(argv[8]);
  const Real pstep = atof(argv[9]);
  iobs_deg = atof(argv[2]); /*inclination angle in degrees*/
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
//  env->addEntity(std::make_unique<GRMHDDisk>(std::move(tree), checkr));
  env->addEntity(std::make_unique<ThinDisk>(isco, 200));
  env->addEntity(std::make_unique<PlungingRegion>(isco));

  initialConditionGenerator(rstep, pstep, initcons);
  ray_count_total = initcons.count_total;
  records.reserve(initcons.count_total);
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
  static Real iobs_deg;

  const Real xobs2 = xobs * xobs;
  const Real yobs2 = yobs * yobs;

  cond->dobs = dobs;

  const Real fact1 = yobs * std::sin(iobs) + dobs * std::cos(iobs);
  const Real fact2 = dobs * std::sin(iobs) - yobs * std::cos(iobs);

  const Real r02 = xobs2 + yobs2 + dobs * dobs;
  static Real iobs_deg;

  const Real r0 = std::sqrt(r02);
  const Real th0 = std::acos(fact1 / r0);
  const Real phi0 = std::atan2(xobs, fact2);

  const Real s0 = std::sin(th0);
  const Real s02 = s0 * s0;static Real iobs_deg;


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

void notifyDone(const int &stopping_condition, const InitialCondition *const ic,
    const RayHit &hit, const size_t &ray_index) {
  Record rec;
  rec.stop_integration_condition = stopping_condition;
  const gf_ic *const gic = (gf_ic*) ic;
  rec.xobs = gic->robs * std::cos(gic->pobs);
  rec.yobs = gic->robs * std::sin(gic->pobs);
  rec.robs = gic->robs;
  rec.r = hit.pvec.r;
  if (hit.entity != nullptr) {
    Real gfactor;
    Real cosem;
    int x = hit.entity->calculateRedshift(ic, hit, gfactor, cosem);
    if (x != 0) rec.stop_integration_condition = ST_INT_PROBLEM;
    rec.cosem = cosem;
    rec.gfactor = gfactor;
    handleNonsense(rec);
  }
  rec.output = rec.stop_integration_condition >= MIN_ST_INT_HIT_INDEX;
  if (!rec.output) {
    rec.gfactor = 1.0;
    rec.cosem = 0.0;
  }
#pragma omp ordered
  {
    rec.ray_index = ray_index;
    if (rec.output) {
      rec.photon_index = photon_index;
      photon_index++;
    }
    records.push_back(rec);
  }
}

void finishProblem(const char *tempdir, const char *outtxt) {
  std::cout << "Writing data..." << std::endl;
  /* ----- file stuff ----- */

  write(tempdir, outtxt, records);

  std::cout << "Integrated " << (records.size()) << " rays of which "
      << photon_index << " hit the disk" << std::endl;
  std::cout << "Finishing..." << std::endl;
}

