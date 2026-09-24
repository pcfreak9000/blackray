#pragma once

#include "def.hpp"

extern Real epsi3, a13, a22, a52;
extern Real spin;

void metric(Real z1, Real z2, Real mn[][4]);
void metric_rderivatives(Real z1, Real z2, Real dmn[][4]);
void scalarProduct(Real met[4][4], Real *fvec0, Real *fvec1, Real &scal);
void correct4VelNorm(Real met[4][4], Real norm, Real *fvel);
