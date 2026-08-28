#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <cmath>

#include "def.hpp"
#include "quadtree.hpp"
#include "raytracingnew.hpp"
#include "find_isco.hpp"

Real epsi3, a13, a22, a52;
Real spin;
Real iobs_deg;
int phicount;

int main(int argc, char *argv[]) {
  std::cout << "Setting up raytracer..." << std::endl;
  if (DEBUG_DIV != 1.0) std::cout << "debug_div is non-one" << std::endl;
  Real spin2;
  Real E_line, N_0, N_tot, alpha;
  Real iobs;
  Real robs;
  Real robs_i, robs_f, rstep, rstep2, pstep;
  Real xin, xout;

  Real isco;
  Real gfactor;

  Real E_obs[IMAX];
  Real N_obs[IMAX];
  Real fphi[IMAX];

  int photon_index = 0;

  char filename_o[256];
  char filename_o2[256];

  FILE *foutput;
  FILE *foutput_coord;

  const char *diskdatafile = argv[10];

  Real maxx;
  Real maxy;
  QuadTree *tree = readFileToTree(diskdatafile, maxx, maxy);
  if (!tree) return 1;

  /* ----- Set free parameters ----- */

  // Input parameters: spin, incl, a13, a22, a52, epsi3, alpha, rstep, pstep
  spin = atof(argv[1]);
  iobs_deg = atof(argv[2]); /*inclination angle in degrees*/
  a13 = atof(argv[3]); /* deformation parameters */
  a22 = atof(argv[4]);
  a52 = atof(argv[5]);
  epsi3 = atof(argv[6]);
  alpha = atof(argv[7]);
  rstep = atof(argv[8]);
  pstep = atof(argv[9]);

  spin2 = spin * spin;

  Real maxr_xdir = std::sqrt(SQR(maxx+10.0) - spin2) + 10.0;
  Real maxr_ydir = std::sqrt(SQR(maxy + 10.0)) + 10.0;
  Real checkr = maxr_ydir;
  if (maxr_xdir > maxr_ydir) checkr = maxr_xdir;
  find_isco(15.0, isco); /* Depends upon the properties of BH */
  GRMHDDisk disk(tree, checkr);
  ThinDisk thin(isco, 100);
  PlungingRegion plunge(isco);
  Env env;
  env.addEntity(&thin);
  env.addEntity(&plunge);
  iobs = Pi / 180 * iobs_deg; /* inclination angle of the observer in rad */
  // iobs = acos(iobs_deg);

  /* ----- Set model for the spectral line ----- */

  E_line = 6.4; /* energy rest of the line in keV */
  N_0 = 1.0; /* normalization */
  // alpha  = -3;     radial power law index

  /* ----- Set inner and outer radius of the disk ----- */

  /*------------------------------------------*/

  /*** thin disk parameters  ***/
  xin = isco; /* inner radius of the accretion disk; set isco */
  xout = 10; /* outer radius of the accretion disk */

  /* ----- Set computational parameters ----- */

  robs_i = 0.1; //this was 1, but that is too big and will yield artifacts
  robs_f = 11;

  // rstep  = 1.008;
  rstep2 = (rstep - 1) / rstep;
  // pstep  = 2*Pi/720;

  E_obs[0] = 0.0125000002; /* minimum photon energy detected by the observer; in keV */
  N_obs[0] = 0;
  Real E_step = 0.025;
  for (size_t i = 1; i < IMAX; i++) {
    E_obs[i] = E_obs[i - 1] + E_step;
    N_obs[i] = 0;
  }
  const char *tempdir = argv[11];
  const char *outtxt = argv[12];
  std::string tempdirs(tempdir);
  /*Iron line output file*/
  // sprintf(filename_o,"iron_a%.03f.epsilon_r%.02f.epsilon_t%.02f.i%.02f.dat",spin,epsi3,iobs_deg);
  // sprintf(filename_o,"ironline_data/iron_a%.05Le.i%.02Le.e_%.02Le.a13_%.02Le.a22_%.02Le.a52_%.02Le.dat",spin,iobs_deg,epsi3,a13,a22,a52);
  std::string s1("ironline_data/"
      "iron_a_%.05Lf_i_%.05Lf_e_%.05Lf_a13_%.05Lf_a22_%.05Lf_a52_%.05Lf.dat");
  snprintf(filename_o, sizeof(filename_o), (tempdirs + s1).c_str(), spin,
      iobs_deg, epsi3, a13, a22, a52);

  /*photon data output file*/
  // sprintf(filename_o2,"coord_a%.03f.epsilon_r%.02f.epsilon_t%.02f.i%.02f.dat",spin,epsi3,iobs_deg);
  std::string s2(
      "data/"
          "photons_data_a%.05Lf_i_%.05Lf_e_%.05Lf_a13_%.05Lf_a22_%.05Lf_a52_%.05Lf.dat");
  snprintf(filename_o2, sizeof(filename_o2), (tempdirs + s2).c_str(), spin,
      iobs_deg, epsi3, a13, a22, a52);
//  cout << tempdirs+s2 << endl;
//  cout << filename_o2 << endl;
  foutput_coord = fopen(filename_o2, "w");
  if (foutput_coord == nullptr) std::cerr << "Problems with data file!"
      << std::endl;
  //string s3("output.txt");
  std::ofstream tmpOutFile(outtxt);
  std::cout << "Starting raytracing loop" << std::endl;
  unsigned long long raycount = 0;
  unsigned long long hitraycount = 0;

  std::vector<Real> pobs_vec;
  for (Real pobs = 0; pobs < 2 * Pi - 0.5 * pstep; pobs = pobs + pstep) {
    pobs_vec.push_back(pobs);
  }
  size_t pobs_count = pobs_vec.size();
  Real *pobs_data = pobs_vec.data();

  std::vector<Real> robs_vec;
  for (robs = robs_i; robs < robs_f; robs = robs * rstep) {
    robs_vec.push_back(robs);
  }
  size_t robs_count = robs_vec.size();
  Real *robs_data = robs_vec.data();

  /* ----- assign photon position in the grid ----- */
  for (size_t robs_index = 0; robs_index < robs_count; robs_index++) {
    Real robs = robs_data[robs_index];
    //for (robs = robs_i; robs < robs_f; robs = robs * rstep) {
    std::cout << "Raytracing: " << (robs - robs_i) / (robs_f - robs_i)
        << std::endl;
    for (size_t i = 0; i < IMAX; i++)
      fphi[i] = 0;

#pragma omp parallel
    {
#pragma omp for ordered schedule(dynamic)
      for (size_t pobs_index = 0; pobs_index < pobs_count; pobs_index++) {
        Real pobs = pobs_data[pobs_index];
        //for (Real pobs = 0; pobs < 2 * Pi - 0.5 * pstep; pobs = pobs + pstep) {
        Real xobs, yobs;
        xobs = robs * std::cos(pobs);
        yobs = robs * std::sin(pobs);
        /*entering in raytrace_new.cpp*/
        // printf("entering in the raytrace part of the code\n");
        int stop_integration_condition = 0;
        RayHit hit;
        raytrace(xobs, yobs, iobs, xin, xout, hit, stop_integration_condition,
            &env);
#pragma omp ordered
        {
          raycount++;
          if ((stop_integration_condition >= 128
              && stop_integration_condition <= 131)
              || stop_integration_condition == 512
              || stop_integration_condition == 600) {
            hitraycount++;
            fprintf(foutput_coord, "%d %Lf %Lf %Lf %Lf %Lf\n", photon_index,
                xobs, yobs, hit.r, hit.gfactor, hit.cosem);

            photon_index++;

            gfactor = hit.gfactor;
            if (!RESTRICT_DEBUGFILE_CRIT && raycount % DEBUGFILE_OUT_DIV == 0) {
              tmpOutFile << xobs << " " << yobs << " " << gfactor << " "
                  << stop_integration_condition << " " << std::endl;
            }
            /* --- integration - part 1 --- */
            Real pp = gfactor * E_line;
            Real bucket = (pp - E_obs[0]) / E_step;
            if (bucket >= 0.0) {
              size_t index = std::floor(bucket);
              if (index < IMAX) {
                Real qq = gfactor * gfactor * gfactor * gfactor;
                qq = qq * std::pow(hit.r, alpha);
                fphi[index] = fphi[index] + qq;

              }
            }
          } else {
            if ((!RESTRICT_DEBUGFILE_CRIT && raycount % DEBUGFILE_OUT_DIV == 0)
                || (stop_integration_condition == 255
                    || stop_integration_condition == 6)) {
              tmpOutFile << xobs << " " << yobs << " " << 1.0 << " "
                  << stop_integration_condition << std::endl;
            }
          }
        }
      }
      /* --- integration - part 2 --- */
#pragma omp for
      for (size_t i = 0; i < IMAX; i++) {
        Real fr = SQR(robs) * fphi[i] * rstep2;
        N_obs[i] += fr;
      }
    }
  }

  N_tot = 0.0;
#pragma omp parallel for reduction(+:N_tot)
  for (size_t i = 0; i < IMAX; i++) {
    N_obs[i] = N_0 * N_obs[i] / E_obs[i];
    N_tot += N_obs[i];
  }

  std::cout << "Integrated " << raycount << " rays of which " << hitraycount
      << " hit the disk" << std::endl;
  std::cout << "Finishing..." << std::endl;
  tmpOutFile.close();
  /* --- print spectrum --- */

  foutput = fopen(filename_o, "w");
  if (foutput == nullptr) std::cerr << "Problems with iron line file!"
      << std::endl;

  for (size_t i = 0; i < IMAX; i++) {
    fprintf(foutput, "%Lf %.10Lf\n", E_obs[i], N_obs[i] / N_tot);
  }

  fclose(foutput);
  fclose(foutput_coord);
  std::cout << "Done" << std::endl;
  return 0;
}
