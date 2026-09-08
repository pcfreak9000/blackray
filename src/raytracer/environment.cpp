#include <cmath>
#include <iostream>

#include "environment.hpp"
#include "metric.hpp"
#include "redshift.hpp"
#include "problem.hpp"





void Env::addEntity(std::unique_ptr<Entity> ptr) {
  maxr = std::max(maxr, ptr->getMaxRadius());
  this->ents.push_back(std::move(ptr));
}

int Env::checkIntersect(const Real &r, const Real &th, const Real &rprev,
    const Real &thprev, Entity *&hitent) {
  if (r > maxr) {
    hitent = nullptr;
    return 0;
  }
  for (std::unique_ptr<Entity> &ent : this->ents) {
    int n = ent->checkIntersect(r, th, rprev, thprev);
    if (n) {
      hitent = ent.get();
      return n;
    }
  }
  hitent = nullptr;
  return 0;
}

int Entity::calculateRedshift(const InitialCondition *ic, const RayHit &hit,
    Real &gfactor, Real &cosem) {
  gfactor = 1.0;
  cosem = 0.0;
  return 0;
}

GRMHDDisk::GRMHDDisk(QuadTree *tree, Real checkr) :
    tree(tree), checkr(checkr) {
}


int GRMHDDisk::checkIntersect(const Real &r, const Real &th, const Real &rprev,
    const Real &thprev) {
  //not at all close to disk so we don't need to perform the checks below
  if (r > checkr) return 0;
  //check if the new position intersects the accretion disk
  //convert coordinates of current and previous position via a BL-cartesian conversion
  Real spin2 = SQR(spin);
  Real xcoord = std::sqrt(r * r + spin2) * std::sin(th);
  Real ycoord = r * std::cos(th);
  Real xcoordprev = std::sqrt(rprev * rprev + spin2) * std::sin(thprev);
  Real ycoordprev = rprev * std::cos(thprev);
  SurfacePoint outSurface;
  bool res = tree->get_interpolated_sp(xcoordprev, ycoordprev, xcoord, ycoord,
      outSurface);
  if (std::abs(outSurface.y) <= 0.0 || outSurface.x <= 0.0) return 0;
  return res ? outSurface.index : 0;
}



int GRMHDDisk::calculateRedshift(const InitialCondition* ic, const RayHit &hit, Real &gfactor, Real &cosem) {
  //to calculate the redshift, we need the photon momentum k (which is present with kr and kth, kt=-E=kt0, kphi=L=kphi0) the observer 4-vel,
  //which is (1,0,0,0), and the interpolated 4-vel of the disk. With this, we can calculate the gfactor.
  Real met[4][4];

  metric(hit.pvec.r, hit.pvec.th, met);
  Real spin2 = SQR(spin);
  Real xcoord = std::sqrt(hit.pvec.r * hit.pvec.r + spin2) * std::sin(hit.pvec.th);
  Real ycoord = hit.pvec.r * std::cos(hit.pvec.th);
  Real xcoordprev = std::sqrt(hit.pvecau.r * hit.pvecau.r + spin2)
      * std::sin(hit.pvecau.th);
  Real ycoordprev = hit.pvecau.r * std::cos(hit.pvecau.th);
  SurfacePoint outSurface;
  tree->get_interpolated_sp(xcoordprev, ycoordprev, xcoord, ycoord, outSurface);
  //Real x = std::std::sqrt(r);
  //Real p_ut = (0.0 + CUBE(x))/std::std::sqrt(CUBE(x)*(2*0.0+CUBE(x)-3*x));
  //Real p_uph = 1/std::std::sqrt(CUBE(x)*(2*0.0+CUBE(x)-3*x));
  //Real uarray[4] = {p_ut,0.0,0.0,p_uph};
  //Real uarray[4] = {1,0,0,0};
//  Real uarray[4] = { outSurface.u0, outSurface.u1, outSurface.u2,
//      outSurface.u3 };

  //to fix any inconsistencies introduced by linear interpolation or the change of coordinate chart (KS->BL)
  //or code differences between Athena++ and Blackray or simply numerical issues in the entire pipeline
  //the fix is done by recalculating (only) the time component of the 4-velocity so that the normalization is correct, i.e. far closer to -1.
  //the highest delta |spi.u0-fixedu0| is approximately 0.05 for an average disk.
  Real norm;
  scalarProduct(met, outSurface.u, outSurface.u, norm);
  correct4VelNorm(met, norm, outSurface.u);
  //uarray[0] is nan sometimes: this should not happen but for some reason the velocity correction introduces this (after code refactorings) so the questions stays, why suddenly now??
  Real newnorm;
  scalarProduct(met, outSurface.u, outSurface.u, newnorm);
#ifdef DEBUG_FVEL_NORM
      if(norm > -0.97 || norm < -1.03) {
        std::cout << "4-Vel norm deviates significantly, ignoring ray" << std::endl;
        std::cout << "Old norm: " << norm << std::endl;
        std::cout << "Fixed norm: " << newnorm << std::endl;
        std::cout << "Delta components: " << spi.u0-uarray[0] << " " << spi.u1-uarray[1] << " " << spi.u2-uarray[2] << " " << spi.u3-uarray[3] << " " << std::endl;
        stop_integration = 6;
      }
#endif
  //if we can't fix this mess, invalidate ray
  if (newnorm > -0.95 || newnorm < -1.05) {
    std::cout << "even fixed 4-vel norm deviates significantly, ignoring ray"
        << std::endl;
    //stop_integration = 6;
    cosem = 0.0;
    gfactor = 1.0;
    return ST_INT_PROBLEM;
  }

  Real g_tt, g_pp, g_tp;
  g_tt = met[0][0];
  g_pp = met[3][3];
  g_tp = met[0][3];
  Real denom = (g_tt * g_pp - g_tp * g_tp);
  Real ktcalc = -(g_pp + ic->b * g_tp) / denom;
  Real kphicalc = (g_tp + ic->b * g_tt) / denom;

  Real karray[4] = { ktcalc, hit.pvec.kr, hit.pvec.kth, kphicalc };

#ifdef DEBUG_FMOM_NORM
        Real knorm;
        scalarProduct(met, karray, karray, knorm);
        if(knorm > 0.03 || knorm < -0.03) {
          std::cout << "4-Momentum norm deviates significantly" << std::endl;
          std::cout << "Norm: " << knorm << std::endl;
          stop_integration = 6;
        }
#endif

  Real emenergy;
  scalarProduct(met, outSurface.u, karray, emenergy);
  gfactor = ((gf_ic*) ic)->obsenergy / emenergy;

  //cosem stays artifical
  Real gfactorforcosem;
  redshift(hit.pvec.r, ((gf_ic*) ic)->const1, gfactorforcosem);
  /*Non Kerr PRD 90, 064002 (2014) Eq. 34*/
  cosem = ((gf_ic*) ic)->carter * gfactorforcosem
      / std::sqrt(SQR(hit.pvec.r) + epsi3 / hit.pvec.r);
  //Workaround for redshift function giving nan...
  if (std::isnan(cosem)) {
    cosem = ((gf_ic*) ic)->carter * gfactor
        / std::sqrt(SQR(hit.pvec.r) + epsi3 / hit.pvec.r);
    if (cosem > 1.05) {
      std::cout << "Cosem was nan, then fixed cosem was > 1.05, ignoring ray: "
          << cosem << std::endl;
      gfactor = 1.0;
      cosem = 0.0;
      return ST_INT_PROBLEM;
    } else if (cosem > 1.0) {
      cosem = 1.0;
    }
  }

  return ST_INT_CONTINUE;
}

Real GRMHDDisk::getMaxRadius() {
  return checkr;
}

ThinDisk::ThinDisk(Real innerr, Real outerr) :
    innerr(innerr), outerr(outerr) {
}

int ThinDisk::checkIntersect(const Real &r, const Real &th, const Real &rprev,
    const Real &thprev) {
  if (r > outerr && rprev > outerr) return false;
  if (r <= innerr && rprev <= innerr) return false;
  if ((th > Pi / 2.0 && thprev < Pi / 2.0)
      || (th < Pi / 2.0 && thprev > Pi / 2.0)) {
    return 512;
  }
  return 0;
}


int ThinDisk::calculateRedshift(const InitialCondition* ic, const RayHit &hit,
    Real &gfactor, Real &cosem) {
  redshift(hit.pvec.r, ((gf_ic*) ic)->const1, gfactor);
  cosem = ((gf_ic*) ic)->carter * gfactor
      / std::sqrt(SQR(hit.pvec.r) + epsi3 / hit.pvec.r);
  return 0;
}

Real ThinDisk::getMaxRadius() {
  return outerr;
}

PlungingRegion::PlungingRegion(Real isco) :
    isco(isco) {
}

int PlungingRegion::checkIntersect(const Real &r, const Real &th,
    const Real &rprev, const Real &thprev) {
  if (r > isco && rprev > isco) return false;
  if ((th > Pi / 2.0 && thprev < Pi / 2.0)
      || (th < Pi / 2.0 && thprev > Pi / 2.0)) {
    return 600;
  }
  return 0;
}
int PlungingRegion::calculateRedshift(const InitialCondition* ic,
    const RayHit &hit, Real &gfactor, Real &cosem) {
  Real met[4][4];
  metric(hit.pvec.r, hit.pvec.th, met);
  Real g_tt, g_pp, g_tp;
  g_tt = met[0][0];
  g_pp = met[3][3];
  g_tp = met[0][3];
  Real denom = (g_tt * g_pp - g_tp * g_tp);
  Real ktcalc = -(g_pp + ic->b * g_tp) / denom;
  Real kphicalc = (g_tp + ic->b * g_tt) / denom;

  Real karray[4] = { ktcalc, hit.pvec.kr, hit.pvec.kth, kphicalc };
  redshift_plunge(isco, hit.pvec.r, karray, gfactor);
  cosem = ((gf_ic*) ic)->carter * gfactor
      / std::sqrt(SQR(hit.pvec.r) + epsi3 / hit.pvec.r);
  return 0;
}

Real PlungingRegion::getMaxRadius() {
  return isco;
}

