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
void getRecord(Record& rec, const RayHit& hit, const InitialCondition*const ic);
void postRecord(Record& rec, const InitialCondition*const ic);
