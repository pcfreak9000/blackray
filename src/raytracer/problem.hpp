#pragma once
#include <memory>
#include "def.hpp"
#include "raytracingnew.hpp"

class gf_ic : public InitialCondition {
public:
  Real carter;
  Real const1;
  Real obsenergy;
  Real robs;
  Real pobs;
};

struct myinitialconditions {
  std::vector<Real> pobs_vec;
  std::vector<Real> robs_vec;
  size_t count_total;
};

void setupProblem(int argc, char *argv[], Env *env,size_t& ray_count_total);
std::unique_ptr<InitialCondition> initialcondition(const size_t& ray_index);
void notifyDone(const int& stopping_condition, const InitialCondition*const ic, const RayHit& hit, const size_t& ray_index);
void finishProblem(const char* tempdir, const char* outtxt);
