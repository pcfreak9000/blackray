#pragma once

#include "def.hpp"
#include "quadtree.hpp"

void raytrace(Real xobs, Real yobs, Real iobs,
    Real rin, Real disk_length_combined, RayHit &hit,
    int &stop_integration, QuadTree* tree, Real checkr);
