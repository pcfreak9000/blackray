#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <cmath>
#include <memory>

#include "def.hpp"
#include "quadtree.hpp"
#include "raytracingnew.hpp"
#include "find_isco.hpp"

Real epsi3, a13, a22, a52;
Real spin;
Real iobs_deg;
int phicount;

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

struct initialconditions {
  std::vector<Real> pobs_vec;
  std::vector<Real> robs_vec;
  size_t count_total;
};

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

void initialConditionGenerator(const Real &rstep, const Real &pstep, initialconditions &ic) {
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

int main(int argc, char *argv[]) {
  std::cout << "Setting up raytracer..." << std::endl;
  if (DEBUG_DIV != 1.0) std::cout << "debug_div is non-one" << std::endl;
  /* ----- Set free parameters ----- */

  // Input parameters: spin, incl, a13, a22, a52, epsi3, alpha, rstep, pstep
  spin = atof(argv[1]);
  iobs_deg = atof(argv[2]); /*inclination angle in degrees*/
  a13 = atof(argv[3]); /* deformation parameters */
  a22 = atof(argv[4]);
  a52 = atof(argv[5]);
  epsi3 = atof(argv[6]);
  const Real alpha = atof(argv[7]);
  const Real rstep = atof(argv[8]);
  const Real pstep = atof(argv[9]);
  const char *diskdatafile = argv[10];
  const char *tempdir = argv[11];
  const char *outtxt = argv[12];

  const Real spin2 = SQR(spin);
  const Real iobs = Pi / 180 * iobs_deg; /* inclination angle of the observer in rad */
  // iobs = acos(iobs_deg);

  /* ----- SETUP ENVIRONMENT ----- */
  Real maxx;
  Real maxy;
  std::unique_ptr<QuadTree> tree = readFileToTree(diskdatafile, maxx, maxy);
  if (!tree) return 1;
  Real maxr_xdir = std::sqrt(SQR(maxx+10.0) - spin2) + 10.0;
  Real maxr_ydir = std::sqrt(SQR(maxy + 10.0)) + 10.0;
  Real checkr = maxr_ydir;
  if (maxr_xdir > maxr_ydir) checkr = maxr_xdir;
  Real isco;
  find_isco(15.0, isco); /* Depends upon the properties of BH */
  GRMHDDisk disk(tree.get(), checkr);
  ThinDisk thin(isco, 200);
  PlungingRegion plunge(isco);
  Env env;
  env.addEntity(&disk);
  //env.addEntity(&thin);
  //env.addEntity(&plunge);



  initialconditions initcons;
  initialConditionGenerator(rstep, pstep, initcons);

  std::vector<Record> records;
  records.reserve(initcons.count_total);

  size_t photon_index = 0;
  std::cout << "Starting raytracing loop" << std::endl;
#pragma omp parallel for ordered schedule(dynamic)
  for (size_t ray_index = 0; ray_index < initcons.count_total; ray_index++) {
    /* ----- unfold initial conditions ----- */
    size_t robs_index = ray_index / initcons.pobs_vec.size();
    size_t pobs_index = ray_index % initcons.pobs_vec.size();
    if (ray_index % 10000 == 0) {
      std::cout << "Progress: " << ray_index / (double) (initcons.count_total)
          << std::endl;
    }
    Real robs = initcons.robs_vec[robs_index];
    Real pobs = initcons.pobs_vec[pobs_index];
    Real xobs = robs * std::cos(pobs);
    Real yobs = robs * std::sin(pobs);

    /* ----- raytrace ----- */
    int stop_integration_condition = 0;
    RayHit hit;
    raytrace(xobs, yobs, iobs, hit, stop_integration_condition, &env);
    bool yes_hit = (stop_integration_condition >= 128
        && stop_integration_condition <= 131)
        || stop_integration_condition == 512
        || stop_integration_condition == 600;
    Record rec;
    rec.stop_integration_condition = stop_integration_condition;
    rec.cosem = hit.cosem;
    rec.gfactor = hit.gfactor;
    rec.r = hit.r;
    rec.xobs = xobs;
    rec.yobs = yobs;
    rec.output = yes_hit;
    rec.robs = robs;
#pragma omp ordered
    {
      rec.ray_index = ray_index;
      if (yes_hit) {
        rec.photon_index = photon_index;
        photon_index++;
      }
      records.push_back(rec);
    }
  }
  std::cout << "Writing data..." << std::endl;
  /* ----- file stuff ----- */

  write(tempdir,outtxt,records);


  std::cout << "Integrated " << (initcons.count_total) << " rays of which "
      << photon_index << " hit the disk" << std::endl;
  std::cout << "Finishing..." << std::endl;

  /* ----- calculate iron line stuff ----- */
  ironlinestuff(rstep, alpha, tempdir, records);


  std::cout << "Done" << std::endl;
  return 0;
}
