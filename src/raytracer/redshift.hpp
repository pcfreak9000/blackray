#pragma once

#include "def.hpp"

void redshift(Real r, Real ktkp, Real &gg);
void redshift_plunge(const Real& isco, const Real& r, Real* kvec, Real& gg);
