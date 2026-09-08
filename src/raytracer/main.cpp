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
#include "problem.hpp"

Real epsi3, a13, a22, a52;
Real spin;
Real iobs_deg;
int phicount;





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
  const char *tempdir = argv[11];
  const char *outtxt = argv[12];

  // iobs = acos(iobs_deg);


  Env env;
  size_t ray_count_total;
  setupProblem(argc,argv,&env, ray_count_total);



  std::vector<Record> records;
  records.reserve(ray_count_total);

  size_t photon_index = 0;
  std::cout << "Starting raytracing loop" << std::endl;
#pragma omp parallel for ordered schedule(dynamic)
  for (size_t ray_index = 0; ray_index < ray_count_total; ray_index++) {
    /* ----- unfold initial conditions ----- */
    if (ray_index % 10000 == 0) {
      std::cout << "Progress: " << ray_index / (double) (ray_count_total)
          << std::endl;
    }



    /* ----- raytrace ----- */
    Record rec;
    rec.gfactor = 1.0;
    rec.cosem = 0.0;
    raytrace(ray_index, rec, &env);
    int& stop_integration_condition = rec.stop_integration_condition;
    bool yes_hit = (stop_integration_condition >= 128
        && stop_integration_condition <= 131)
        || stop_integration_condition == 512
        || stop_integration_condition == 600;
    rec.output = yes_hit;
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


  std::cout << "Integrated " << (ray_count_total) << " rays of which "
      << photon_index << " hit the disk" << std::endl;
  std::cout << "Finishing..." << std::endl;

  /* ----- calculate iron line stuff ----- */
  ironlinestuff(rstep, alpha, tempdir, records);


  std::cout << "Done" << std::endl;
  return 0;
}
