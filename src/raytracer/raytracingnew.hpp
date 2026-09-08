#pragma once

#define ST_INT_PHOTON_HORIZON_CROSS_DCRIT 4
#define ST_INT_PHOTON_HORIZON_CROSS_RCRIT 5
#define ST_INT_PROBLEM 6
#define ST_INT_PHOTON_ESCAPES 7
#define ST_INT_MAX_ITERATIONS 255
#define ST_INT_CONTINUE 0

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
};

struct RayHit {
  PhaseVec pvec;
  PhaseVec pvecau;
  Entity *entity;
};

void raytrace(const size_t& ray_index, Record& rec, Env *env);
