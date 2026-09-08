#pragma once

#define ST_INT_PHOTON_HORIZON_CROSS_DCRIT 4
#define ST_INT_PHOTON_HORIZON_CROSS_RCRIT 5
#define ST_INT_PROBLEM 6
#define ST_INT_PHOTON_ESCAPES 7
#define ST_INT_MAX_ITERATIONS 1
#define ST_INT_CONTINUE 0

#define MIN_ST_INT_HIT_INDEX 128

class InitialCondition;
struct RayHit;

#include "def.hpp"
#include "environment.hpp"

union PhaseVec {
  Real vars[5];
  struct {
    Real r, th, phi;
    Real kr, kth;
  };
};

class InitialCondition {
public:
  PhaseVec pvec0;
  Real b;
  Real dobs;
};

struct RayHit {
  PhaseVec pvec;
  PhaseVec pvecau;
  Entity *entity;
};
void raytrace(Env *env, size_t ray_count_total);
void raytrace(const size_t& ray_index, Env *env);
