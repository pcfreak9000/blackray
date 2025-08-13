#include "def.hpp"
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>



QuadTree::QuadTree(Real x, Real y, Real width, Real height) :
    x(x), y(y), width(width), height(height), max_elements(4), is_leaf(true) {
}

Real QuadTree::check_intersect(Real x1, Real y1, Real x2, Real y2, SurfaceElement*& out){
  for(SurfaceElement* elem : myelements) {
    Real res = checkIntersect(x1,y1,x2,y2,elem->sp0->x,elem->sp0->y,elem->sp1->x,elem->sp1->y);
    if(res != NO_INTERSECT) {
      out = elem;
      return res;
    }
  }
  if(!is_leaf) {
    for(int i=0; i<4; i++){
      Real ressub = subtrees[i]->check_intersect(x1,y1,x2,y2,out);
      if(ressub != NO_INTERSECT){
        return ressub;
      }
    }
  }
  return NO_INTERSECT;
}

bool QuadTree::fits(SurfaceElement *element) {
  Real maxx = std::max(element->sp0->x, element->sp1->x);
  Real minx = std::min(element->sp0->x, element->sp1->x);
  Real maxy = std::max(element->sp0->y, element->sp1->y);
  Real miny = std::min(element->sp0->y, element->sp1->y);

  return maxx < x+width && minx >= x && maxy < y+height && miny >= y;
}

void QuadTree::subdivide() {
  QuadTree* q1 = new QuadTree(x+width/2.0,y+height/2.0,width/2.0,height/2.0);
  QuadTree* q2 = new QuadTree(x,y+height/2.0,width/2.0,height/2.0);
  QuadTree* q3 = new QuadTree(x,y,width/2.0,height/2.0);
  QuadTree* q4 = new QuadTree(x+width/2.0,y,width/2.0, height/2.0);
  is_leaf = false;
  subtrees.push_back(q1);
  subtrees.push_back(q2);
  subtrees.push_back(q3);
  subtrees.push_back(q4);
  vector<size_t> forremoval;
  for(size_t k=0; k<myelements.size(); k++) {
    SurfaceElement *se = myelements[k];
    for(int i=0; i<4; i++){
      if(subtrees[i]->fits(se)){
        forremoval.push_back(k);
        subtrees[i]->put_element(se);
        break;
      }
    }
  }
  for(size_t rem : forremoval) {
    myelements.erase(myelements.begin()+rem);
  }
}

void QuadTree::put_element(SurfaceElement *element){
  if(is_leaf) {
    myelements.push_back(element);
    if(myelements.size() > max_elements){
      subdivide();
    }
  } else {
    for(int i=0; i<4; i++){
      if(subtrees[i]->fits(element)){
        subtrees[i]->put_element(element);
        break;
      }
      if(i==3){
        myelements.push_back(element);
      }
    }
  }
}

int main(int argc, char *argv[]) {
  std::cout << "Setting up raytracer..." << std::endl;
  long double spin2;
  long double E_line, N_0, N_tot, N_tot1, N_tot2, alpha;
  long double iobs, dobs;
  long double robs, pobs;
  long double robs_i, robs_f, rstep, rstep2, pstep;
  long double xin, xout;
  long double pp, qq;
  long double fr;

  long double isco;
  long double gfactor;

  long double E_obs[imax];
  long double N_obs[imax];
  long double fphi[imax], fphi0[imax];
  RayHit hit;

  int stop_integration_condition = 0;
  int n1, n2, n3;
  int i, j, m;
  int photon_index = 0;

  char filename_i[128];
  char filename_o[128];
  char filename_o2[128];

  FILE *finput;
  FILE *foutput;
  FILE *foutput_coord;

  const char *diskdatafile = argv[10];
  std::ifstream disk(diskdatafile);
  if (!disk) {
    std::cerr << "Error: Could not open file!" << std::endl;
    return 1;
  }

  std::string line;
  std::getline(disk, line); // Read and ignore the header line

  std::vector<SurfacePoint> diskdata, dd2;
  Real minx=1e9,maxx=0,miny=1e9,maxy=0;
  while (std::getline(disk, line)) {
    std::istringstream iss(line);
    SurfacePoint dp,dpunder;

    if (!(iss >> dp.x >> dp.y >> dp.density >> dp.u0 >> dp.u1 >> dp.u2 >> dp.u3)) {
      std::cerr << "Error: Malformed line - " << line << std::endl;
      continue;
    }
    dp.y = dp.y < 0.0 ? 0.0 : dp.y;
    dpunder = dp;
    dpunder.y *= -1;
    if(dp.x < minx) minx = dp.x;
    if(dp.y < miny) miny = dp.y;
    if(dp.x > maxx) maxx = dp.x;
    if(dp.y > maxy) maxy = dp.y;
    if(dpunder.y < miny) miny = dpunder.y;
    if(dpunder.y > maxy) maxy = dpunder.y;
    diskdata.push_back(dp);
    dd2.push_back(dpunder);
  }

  disk.close();
  QuadTree tree(minx-1,miny-1, maxx-minx+2, maxy-miny+2);
  for(int i=0; i<diskdata.size()-1; i++){
    SurfaceElement elem;
    elem.sp0 = &diskdata[i];
    elem.sp1 = &diskdata[i+1];
    tree.put_element(&elem);
  }
  for(int i=0; i<dd2.size()-1; i++){
    SurfaceElement elem;
    elem.sp0 = &dd2[i];
    elem.sp1 = &dd2[i+1];
    tree.put_element(&elem);
  }

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

  iobs = Pi / 180 * iobs_deg; /* inclination angle of the observer in rad */
  // iobs = acos(iobs_deg);

  /* ----- Set model for the spectral line ----- */

  E_line = 6.4; /* energy rest of the line in keV */
  N_0 = 1.0; /* normalization */
  // alpha  = -3;     radial power law index

  /* ----- Set inner and outer radius of the disk ----- */

  find_isco(15.0, isco); /* Depends upon the properties of BH */

  /*------------------------------------------*/

  /*** thin disk parameters  ***/
  xin = isco; /* inner radius of the accretion disk; set isco */
  xout = 50; /* outer radius of the accretion disk */

  /* ----- Set computational parameters ----- */

  robs_i = 1;
  robs_f = 80;

  // rstep  = 1.008;
  rstep2 = (rstep - 1) / rstep;
  // pstep  = 2*Pi/720;

  E_obs[0] = 0.0125000002; /* minimum photon energy detected by the observer; in keV */
  N_obs[0] = 0;
  for (i = 1; i <= imax - 1; i++) {
    E_obs[i] = E_obs[i - 1] + 0.025;
    N_obs[i] = 0;
  }

  /*Iron line output file*/
  // sprintf(filename_o,"iron_a%.03f.epsilon_r%.02f.epsilon_t%.02f.i%.02f.dat",spin,epsi3,iobs_deg);
  // sprintf(filename_o,"ironline_data/iron_a%.05Le.i%.02Le.e_%.02Le.a13_%.02Le.a22_%.02Le.a52_%.02Le.dat",spin,iobs_deg,epsi3,a13,a22,a52);
  snprintf(filename_o, sizeof(filename_o), "ironline_data/"
      "iron_a_%.05Lf_i_%.05Lf_e_%.05Lf_a13_%.05Lf_a22_%.05Lf_a52_%.05Lf.dat",
      spin, iobs_deg, epsi3, a13, a22, a52);

  /*photon data output file*/
  // sprintf(filename_o2,"coord_a%.03f.epsilon_r%.02f.epsilon_t%.02f.i%.02f.dat",spin,epsi3,iobs_deg);
  snprintf(filename_o2, sizeof(filename_o2), "data/"
      "photons_data_a%.05Lf_i_%.05Lf_e_%.05Lf_a13_%.05Lf_a22_%.05Lf_a52_%."
      "05Lf.dat", spin, iobs_deg, epsi3, a13, a22, a52);

  foutput_coord = fopen(filename_o2, "w");
  std::ofstream tmpOutFile("output.txt");
  std::cout << "Starting raytracing loop" << std::endl;
  /* ----- assign photon position in the grid ----- */
  for (robs = robs_i; robs < robs_f; robs = robs * rstep) {
    std::cout << "Raytracing: " << (robs - robs_i) / (robs_f - robs_i)
        << std::endl;

    for (i = 0; i <= imax - 1; i++)
      fphi[i] = 0;

    for (pobs = 0; pobs < 2 * Pi - 0.5 * pstep; pobs = pobs + pstep) {
      xobs = robs * cos(pobs);
      yobs = robs * sin(pobs);

      /*entering in raytrace_new.cpp*/

      // printf("entering in the raytrace part of the code\n");
      raytrace(xobs, yobs, iobs, xin, xout, hit, stop_integration_condition,
          diskdata.data(), diskdata.size());

      if (stop_integration_condition == 1
          || stop_integration_condition == 128) {
        fprintf(foutput_coord, "%d %Lf %Lf %Lf %Lf %Lf\n", photon_index, xobs,
            yobs, hit.r, hit.gfactor, hit.cosem);

        photon_index++;

        gfactor = hit.gfactor;
        pp = gfactor * E_line;
        tmpOutFile << xobs << " " << yobs << " " << gfactor << " "
            << stop_integration_condition << " " << hit.hc << std::endl;

        /* --- integration - part 1 --- */

        for (i = 0; i <= imax - 2; i++) {
          if (E_obs[i] < pp && E_obs[i + 1] > pp) {
            qq = gfactor * gfactor * gfactor * gfactor;
            qq = qq * pow(hit.r, alpha);

            fphi[i] = fphi[i] + qq;
          }
        }
      } else {
        tmpOutFile << xobs << " " << yobs << " " << 1.0 << " "
            << stop_integration_condition << " " << hit.hc << std::endl;
      }
    }

    /* --- integration - part 2 --- */

    for (i = 0; i <= imax - 1; i++) {
      fr = robs * robs * fphi[i] * rstep2;
      N_obs[i] = N_obs[i] + fr;
    }
  }
  std::cout << "Finishing..." << std::endl;
  tmpOutFile.close();
  /* --- print spectrum --- */

  foutput = fopen(filename_o, "w");
  N_tot = 0.0;

  for (i = 0; i <= imax - 1; i++) {
    N_obs[i] = N_0 * N_obs[i] / E_obs[i];
    N_tot = N_tot + N_obs[i];
  }

  for (i = 0; i <= imax - 1; i++) {
    fprintf(foutput, "%Lf %.10Lf\n", E_obs[i], N_obs[i] / N_tot);
  }

  fclose(foutput);
  fclose(foutput_coord);
  std::cout << "Done" << std::endl;
  return 0;
}
