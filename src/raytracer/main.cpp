#include "def.hpp"
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>



QuadTree::QuadTree(Real x, Real y, Real width, Real height) :
    x(x), y(y), width(width), height(height), max_elements(4), is_leaf(true),level(0) {
}

Real QuadTree::check_intersect(Real x1, Real y1, Real x2, Real y2, SurfaceElement** out){

  if(!is_leaf) {
    for(int i=0; i<4; i++){
      if(subtrees[i]->overlaps(x1,y1,x2,y2)){
        Real ressub = subtrees[i]->check_intersect(x1,y1,x2,y2,out);
        if(ressub != NO_INTERSECT){
          return ressub;
        }
        //this needs some benchmarking to check if this actually improves things:
        if(subtrees[i]->fully_inside(x1,y1,x2,y2)){
          break;
        }
      }
    }
  }
  for(int i=0; i<myelements.size(); i++) {
        SurfaceElement* elem = myelements[i];
        if (elem->sp0->y==0.0 && elem->sp1->y==0.0){
          //if this raytracer is generalized, this does not really belong here. But we ignore segments with y:=0 of the accretion disk.
          continue;
        }
        Real res = checkIntersect(elem->sp0->x,elem->sp0->y,elem->sp1->x,elem->sp1->y,x1,y1,x2,y2);
        if(res != NO_INTERSECT) {
          *out = elem;
          return res;
        }
    }

  return NO_INTERSECT;
}

size_t QuadTree::size(){
  if(is_leaf) return myelements.size();
  size_t s = myelements.size();
  for(int i=0; i<4; i++){
    s += subtrees[i]->size();
  }
  return s;
}

void QuadTree::validate(){
  //std::cout << "Validating..." << std::endl;
  for(SurfaceElement* se : myelements){
    if(!se) std::cout << level << " " << myelements.size() << "--- uh oh ---: " << se << std::endl;
    if(!se->sp0) std::cout << level << " " << myelements.size() << " oh no 0: " << se->sp0 << std::endl;
    if(!se->sp1) std::cout << level << " " << myelements.size() << " oh no 1: " << se->sp1 << std::endl;
  }
#ifndef abc
  if(!is_leaf){
    for(QuadTree* qt : subtrees){
      qt->validate();
    }
  }
#endif
}

bool QuadTree::fits(SurfaceElement *element) {
  Real maxx = std::max(element->sp0->x, element->sp1->x);
  Real minx = std::min(element->sp0->x, element->sp1->x);
  Real maxy = std::max(element->sp0->y, element->sp1->y);
  Real miny = std::min(element->sp0->y, element->sp1->y);

  return maxx < x+width && minx >= x && maxy < y+height && miny >= y;
}

bool QuadTree::overlaps(Real x0, Real y0, Real x1, Real y1) {
  Real maxx = std::max(x0, x1);
  Real minx = std::min(x0, x1);
  Real maxy = std::max(y0, y1);
  Real miny = std::min(y0, y1);
  return maxx >= x && minx < x+width && maxy >= y && miny < y+height;
}

bool QuadTree::fully_inside(Real x0, Real y0, Real x1, Real y1) {
  Real maxx = std::max(x0, x1);
  Real minx = std::min(x0, x1);
  Real maxy = std::max(y0, y1);
  Real miny = std::min(y0, y1);

  return maxx < x+width && minx >= x && maxy < y+height && miny >= y;
}
void QuadTree::subdivide() {
  Real w_2 = width/2.0;
  Real h_2 = height/2.0;
  QuadTree* q1 = new QuadTree(x+w_2, y+h_2, w_2, h_2);
  QuadTree* q2 = new QuadTree(x,     y+h_2, w_2, h_2);
  QuadTree* q3 = new QuadTree(x,     y,     w_2, h_2);
  QuadTree* q4 = new QuadTree(x+w_2, y,     w_2, h_2);
  q1->level = this->level+1;
  q2->level = this->level+1;
  q3->level = this->level+1;
  q4->level = this->level+1;
  is_leaf = false;
  subtrees.push_back(q1);
  subtrees.push_back(q2);
  subtrees.push_back(q3);
  subtrees.push_back(q4);
  for(auto vecit = myelements.begin(); vecit != myelements.end(); vecit++){
    SurfaceElement* se = *vecit;
    for(int i=0; i<4; i++){
      if(subtrees[i]->fits(se)){
        subtrees[i]->put_element(se);
        myelements.erase(vecit);
        vecit--;
        break;
      }
    }
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
  if(DEBUG_DIV != 1.0) std::cout << "debug_div is non-one" << std::endl;
  long double spin2;
  long double E_line, N_0, N_tot, N_tot1, N_tot2, alpha;
  long double iobs, dobs;
  long double robs;//, pobs;
  long double robs_i, robs_f, rstep, rstep2, pstep;
  long double xin, xout;
  long double pp, qq;
  long double fr;

  long double isco;
  long double gfactor;

  long double E_obs[imax];
  long double N_obs[imax];
  long double fphi[imax];
  long double fphi0[imax];



  int n1, n2, n3;
  int i, j, m;
  int photon_index = 0;

  char filename_i[256];
  char filename_o[256];
  char filename_o2[256];

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

  std::vector<SurfacePoint*> diskdata, dd2;
  Real minx=INFINITY,maxx=-INFINITY,miny=INFINITY,maxy=-INFINITY;
  while (std::getline(disk, line)) {
    std::istringstream iss(line);
    SurfacePoint* dpp = new SurfacePoint;
    SurfacePoint* dpunderp = new SurfacePoint;
    SurfacePoint& dp = *dpp;
    SurfacePoint& dpunder = *dpunderp;
    if (!(iss >> dp.x >> dp.y >> dp.density >> dp.u0 >> dp.u1 >> dp.u2 >> dp.u3)) {
      std::cerr << "Error: Malformed line - " << line << std::endl;
      continue;
    }
    if (dp.y < 0.0) {
      std::cerr << "Error: Negative disk heights not supported, results might be significantly wrong!" << std::endl;
      //dp.y *= -1;
      continue;
    }
    dp.y /= DEBUG_DIV;
    dpunder = dp;
    dpunder.y *= -1;
    if(dp.x < minx) minx = dp.x;
    if(dp.y < miny) miny = dp.y;
    if(dp.x > maxx) maxx = dp.x;
    if(dp.y > maxy) maxy = dp.y;
    if(dpunder.y < miny) miny = dpunder.y;
    if(dpunder.y > maxy) maxy = dpunder.y;
    diskdata.push_back(dpp);
    dd2.push_back(dpunderp);
  }

  disk.close();
  QuadTree* treep = new QuadTree(minx-1,miny-1, maxx-minx+2, maxy-miny+2);
  QuadTree& tree = *treep;
  for(int i=0; i<diskdata.size()-1; i++){
    SurfaceElement* elem = new SurfaceElement;
    elem->sp0 = diskdata[i];
    elem->sp1 = diskdata[i+1];
    elem->index = 128+(i/RING_DIV)%2;
    tree.put_element(elem);
  }
  for(int i=0; i<dd2.size()-1; i++){
    SurfaceElement* elem = new SurfaceElement;
    elem->sp0 = dd2[i];
    elem->sp1 = dd2[i+1];
    elem->index = 130+(i/RING_DIV)%2;
    tree.put_element(elem);
  }
  tree.validate();


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

  Real maxr_xdir = sqrt(SQR(maxx+10.0)-spin2)+10.0;
  Real maxr_ydir = sqrt(SQR(maxy+10.0))+10.0;
  Real checkr = maxr_ydir;
  if(maxr_xdir > maxr_ydir) checkr = maxr_xdir;

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
  xout = 200; /* outer radius of the accretion disk */

  /* ----- Set computational parameters ----- */

  robs_i = 1;
  robs_f = 215;

  // rstep  = 1.008;
  rstep2 = (rstep - 1) / rstep;
  // pstep  = 2*Pi/720;

  E_obs[0] = 0.0125000002; /* minimum photon energy detected by the observer; in keV */
  N_obs[0] = 0;
  for (i = 1; i <= imax - 1; i++) {
    E_obs[i] = E_obs[i - 1] + 0.025;
    N_obs[i] = 0;
  }
  const char *tempdir = argv[11];
  const char *outtxt = argv[12];
  string tempdirs(tempdir);
  /*Iron line output file*/
  // sprintf(filename_o,"iron_a%.03f.epsilon_r%.02f.epsilon_t%.02f.i%.02f.dat",spin,epsi3,iobs_deg);
  // sprintf(filename_o,"ironline_data/iron_a%.05Le.i%.02Le.e_%.02Le.a13_%.02Le.a22_%.02Le.a52_%.02Le.dat",spin,iobs_deg,epsi3,a13,a22,a52);
  string s1("ironline_data/"
      "iron_a_%.05Lf_i_%.05Lf_e_%.05Lf_a13_%.05Lf_a22_%.05Lf_a52_%.05Lf.dat");
  snprintf(filename_o, sizeof(filename_o), (tempdirs+s1).c_str(),
      spin, iobs_deg, epsi3, a13, a22, a52);

  /*photon data output file*/
  // sprintf(filename_o2,"coord_a%.03f.epsilon_r%.02f.epsilon_t%.02f.i%.02f.dat",spin,epsi3,iobs_deg);
  string s2("data/"
      "photons_data_a%.05Lf_i_%.05Lf_e_%.05Lf_a13_%.05Lf_a22_%.05Lf_a52_%.05Lf.dat");
  snprintf(filename_o2, sizeof(filename_o2), (tempdirs+s2).c_str(), spin, iobs_deg, epsi3, a13, a22, a52);
//  cout << tempdirs+s2 << endl;
//  cout << filename_o2 << endl;
  foutput_coord = fopen(filename_o2, "w");
  if(foutput_coord==nullptr) std::cerr << "Problems with data file!" << std::endl;
  //string s3("output.txt");
  std::ofstream tmpOutFile(outtxt);
  std::cout << "Starting raytracing loop" << std::endl;
  unsigned long long raycount = 0;
  unsigned long long hitraycount = 0;
  
  std::vector<long double> pobs_vec;
  for (long double pobs = 0; pobs < 2 * Pi - 0.5 * pstep; pobs = pobs + pstep) {
    pobs_vec.push_back(pobs);
  }
  size_t pobs_count = pobs_vec.size();
  long double* pobs_data = pobs_vec.data();
  
  std::vector<long double> robs_vec;
  for (robs = robs_i; robs < robs_f; robs = robs * rstep) {
    robs_vec.push_back(robs);
  }
  size_t robs_count = robs_vec.size();
  long double* robs_data = robs_vec.data();
  
  /* ----- assign photon position in the grid ----- */
  for(size_t robs_index = 0; robs_index < robs_count; robs_index++){
  long double robs = robs_data[robs_index];
  //for (robs = robs_i; robs < robs_f; robs = robs * rstep) {
    std::cout << "Raytracing: " << (robs - robs_i) / (robs_f - robs_i)
        << std::endl;
    for (i = 0; i <= imax - 1; i++)
      fphi[i] = 0;

    #pragma omp parallel for ordered schedule(dynamic)
    for(size_t pobs_index = 0; pobs_index < pobs_count; pobs_index++){
      long double pobs = pobs_data[pobs_index];
    //for (long double pobs = 0; pobs < 2 * Pi - 0.5 * pstep; pobs = pobs + pstep) {
      long double xobs, yobs;
      xobs = robs * cos(pobs);
      yobs = robs * sin(pobs);
      /*entering in raytrace_new.cpp*/
      // printf("entering in the raytrace part of the code\n");
        int stop_integration_condition = 0;
          RayHit hit;
      raytrace(xobs, yobs, iobs, xin, xout, hit, stop_integration_condition,
          diskdata.data(), diskdata.size(), treep, checkr);
      #pragma omp ordered
      { 
      raycount++;
      if (stop_integration_condition >= 128
          && stop_integration_condition <= 131) {
        hitraycount++;
        fprintf(foutput_coord, "%d %Lf %Lf %Lf %Lf %Lf\n", photon_index, xobs,
            yobs, hit.r, hit.gfactor, hit.cosem);

        photon_index++;

        gfactor = hit.gfactor;
        pp = gfactor * E_line;
        if(!RESTRICT_DEBUGFILE_CRIT && raycount % DEBUGFILE_OUT_DIV == 0){
          tmpOutFile << xobs << " " << yobs << " " << gfactor << " "
              << stop_integration_condition << " " << hit.hc << std::endl;
        }
        /* --- integration - part 1 --- */

        for (i = 0; i <= imax - 2; i++) {
          if (E_obs[i] < pp && E_obs[i + 1] > pp) {
            qq = gfactor * gfactor * gfactor * gfactor;
            qq = qq * pow(hit.r, alpha);

            fphi[i] = fphi[i] + qq;
          }
        }
      } else {
        if((!RESTRICT_DEBUGFILE_CRIT && raycount % DEBUGFILE_OUT_DIV == 0) || ( stop_integration_condition==255 || stop_integration_condition==6 )){
        tmpOutFile << xobs << " " << yobs << " " << 1.0 << " "
            << stop_integration_condition << " " << hit.hc << std::endl;
        }
      }
      }
    }
    /* --- integration - part 2 --- */

    for (i = 0; i <= imax - 1; i++) {
      fr = robs * robs * fphi[i] * rstep2;
      N_obs[i] = N_obs[i] + fr;
    }
  }
//  tree.validate();
  std::cout << "Integrated " << raycount << " rays of which " << hitraycount << " hit the disk" << std::endl;
  std::cout << "Finishing..." << std::endl;
  tmpOutFile.close();
  /* --- print spectrum --- */

  foutput = fopen(filename_o, "w");
  if(foutput==nullptr) std::cerr << "Problems with iron line file!" << std::endl;

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
  tree.validate();
  std::cout << "Done" << std::endl;
  return 0;
}
