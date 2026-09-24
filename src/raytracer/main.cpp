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
#include "metric.hpp"

int main(int argc, char *argv[]) {
  std::cout << "Setting up raytracer..." << std::endl;
  if (DEBUG_DIV != 1.0) std::cout << "debug_div is non-one" << std::endl;
  /* ----- Set free parameters ----- */

  // Input parameters: spin, incl, a13, a22, a52, epsi3, alpha, rstep, pstep
  spin = atof(argv[1]);
  a13 = atof(argv[3]); /* deformation parameters */
  a22 = atof(argv[4]);
  a52 = atof(argv[5]);
  epsi3 = atof(argv[6]);
//  const Real alpha = atof(argv[7]);
//  const Real rstep = atof(argv[8]);
//  const Real pstep = atof(argv[9]);
  const char *tempdir = argv[11];
  const char *outtxt = argv[12];

  // iobs = acos(iobs_deg);


  Env env;
  size_t ray_count_total;
  setupProblem(argc,argv,&env, ray_count_total);

  raytrace(&env, ray_count_total);

  finishProblem(tempdir,outtxt);

  /* ----- calculate iron line stuff ----- */
 // ironlinestuff(rstep, alpha, tempdir, records);


  std::cout << "Done" << std::endl;
  return 0;
}
