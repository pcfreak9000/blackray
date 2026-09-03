#pragma once

#include "def.hpp"
#include "environment.hpp"

void raytrace(Real xobs, Real yobs, Real iobs, RayHit &hit, int &stop_integration, Env* env);
