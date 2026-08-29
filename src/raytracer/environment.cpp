#include "environment.hpp"
#include "metric.hpp"
#include "redshift.hpp"
#include <cmath>
#include <iostream>

Real interpolate(Real a, Real b, Real f) {
  return (1.0 - f) * a + f * b;
}
//finds any intersection. it is not guranteed (but very likely) that it is the closest one, i.e. the first hit. Sufficiently small stepsize should circuvent the problem,
//as well as retaking the step with a smaller stepsize.
//this might interpolate between a mesh cell right beyond the horizon, which might contain invalid data, and the first cell outside the horizon,
//for disk shapes which come very close to the horizon
//at the moment, this does not seem to be a problem, but keep this in mind, especially if and when doing a clean rewrite of this
bool get_interpolated_sp(const Real x1, const Real y1, const Real x2,
    const Real y2, QuadTree *quadtree, SurfacePoint &out) {
  SurfaceElement *elem;
  Real result = quadtree->check_intersect(x1, y1, x2, y2, &elem);
  if (result != NO_INTERSECT) {
    out.index = elem->index;
    Real xi = (elem->sp0->x) + result * ((elem->sp1->x) - (elem->sp0->x));
    Real yi = (elem->sp0->y) + result * ((elem->sp1->y) - (elem->sp0->y));
    out.x = xi;
    out.y = yi;
    out.density = interpolate(elem->sp0->density, elem->sp1->density, result);
    //should be zero if either p is zero because linear interpolation of these velocities at that place is probably not physically
    out.u0 = interpolate(elem->sp0->u0, elem->sp1->u0, result);
    out.u1 = interpolate(elem->sp0->u1, elem->sp1->u1, result);
    out.u2 = interpolate(elem->sp0->u2, elem->sp1->u2, result);
    out.u3 = interpolate(elem->sp0->u3, elem->sp1->u3, result);
    return true;
  }
  return false;
}

GRMHDDisk::GRMHDDisk(QuadTree *tree, Real checkr) :
    tree(tree), checkr(checkr) {
}

bool GRMHDDisk::checkIntersect(const Real &r, const Real &th, const Real &rprev,
    const Real &thprev, SurfacePoint &outSurface) {
  //not at all close to disk so we don't need to perform the checks below
  if (r > checkr) return false;
  //check if the new position intersects the accretion disk
  //convert coordinates of current and previous position via a BL-cartesian conversion
  Real spin2 = SQR(spin);
  Real xcoord = std::sqrt(r * r + spin2) * std::sin(th);
  Real ycoord = r * std::cos(th);
  Real xcoordprev = std::sqrt(rprev * rprev + spin2) * std::sin(thprev);
  Real ycoordprev = rprev * std::cos(thprev);

  bool res = get_interpolated_sp(xcoordprev, ycoordprev, xcoord, ycoord, tree,
      outSurface);
  if (std::abs(outSurface.y) <= 0.0 || outSurface.x <= 0.0) return false;
  return res;
}

void scalarProduct(Real met[4][4], Real *fvec0, Real *fvec1, Real &scal) {
  scal = 0.0;
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      scal += met[i][j] * fvec0[i] * fvec1[j];
    }
  }
}

void correct4VelNorm(Real met[4][4], Real norm, Real *fvel) {
  Real g_tt = met[0][0];
  Real g_tp = met[0][3];
  Real udt = fvel[0];
  Real up = fvel[3];
  Real Fp1 = norm + 1;

  Real radi = SQR(g_tt*udt) + 2.0 * g_tt * udt * g_tp * up
      - g_tt * Fp1+SQR(g_tp*up);
  fvel[0] = -(std::sqrt(radi) + g_tp * up) / g_tt;

  //Real dif = norm+1;
  //Real deltaut = dif/met[0][0];
  //fvel[0] = std::sqrt(fvel[0]*fvel[0]-deltaut);
}
int GRMHDDisk::intersect(const IntegratorData &id,
    const SurfacePoint &surfacepoint, RayHit &hit) {
  hit.r = id.r;
  //to calculate the redshift, we need the photon momentum k (which is present with kr and kth, kt=-E=kt0, kphi=L=kphi0) the observer 4-vel,
  //which is (1,0,0,0), and the interpolated 4-vel of the disk. With this, we can calculate the gfactor.
  Real met[4][4];
  metric(id.r, id.th, met);

  //Real x = std::std::sqrt(r);
  //Real p_ut = (0.0 + CUBE(x))/std::std::sqrt(CUBE(x)*(2*0.0+CUBE(x)-3*x));
  //Real p_uph = 1/std::std::sqrt(CUBE(x)*(2*0.0+CUBE(x)-3*x));
  //Real uarray[4] = {p_ut,0.0,0.0,p_uph};
  //Real uarray[4] = {1,0,0,0};
  Real uarray[4] = { surfacepoint.u0, surfacepoint.u1, surfacepoint.u2,
      surfacepoint.u3 };

  //to fix any inconsistencies introduced by linear interpolation or the change of coordinate chart (KS->BL)
  //or code differences between Athena++ and Blackray or simply numerical issues in the entire pipeline
  //the fix is done by recalculating (only) the time component of the 4-velocity so that the normalization is correct, i.e. far closer to -1.
  //the highest delta |spi.u0-fixedu0| is approximately 0.05 for an average disk.
  Real norm;
  scalarProduct(met, uarray, uarray, norm);
  correct4VelNorm(met, norm, uarray);
  //uarray[0] is nan sometimes: this should not happen but for some reason the velocity correction introduces this (after code refactorings) so the questions stays, why suddenly now??
  Real newnorm;
  scalarProduct(met, uarray, uarray, newnorm);
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
    hit.cosem = 0.0;
    hit.gfactor = 1.0;
    hit.r = id.r;
    return 6;
  } else {

    Real g_tt, g_pp, g_tp;
    g_tt = met[0][0];
    g_pp = met[3][3];
    g_tp = met[0][3];
    Real denom = (g_tt * g_pp - g_tp * g_tp);
    Real ktcalc = -(g_pp + id.b * g_tp) / denom;
    Real kphicalc = (g_tp + id.b * g_tt) / denom;

    Real karray[4] = { ktcalc, id.kr, id.kth, kphicalc };

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
    scalarProduct(met, uarray, karray, emenergy);
    hit.gfactor = id.obsenergy / emenergy;

    //cosem stays artifical
    Real gfactorforcosem;
    redshift(hit.r, id.const1, gfactorforcosem);
    /*Non Kerr PRD 90, 064002 (2014) Eq. 34*/
    hit.cosem = id.carter * gfactorforcosem
        / std::sqrt(SQR(hit.r) + epsi3 / hit.r);
    //Workaround for redshift function giving nan...
    if (std::isnan(hit.cosem)) {
      hit.cosem = id.carter * hit.gfactor
          / std::sqrt(SQR(hit.r) + epsi3 / hit.r);
      if (hit.cosem > 1.05) {
        std::cout
            << "Cosem was nan, then fixed cosem was > 1.05, ignoring ray: "
            << hit.cosem << std::endl;
        hit.gfactor = 1.0;
        hit.cosem = 0.0;
        return 6;
      } else if (hit.cosem > 1.0) {
        hit.cosem = 1.0;
      }
    }
  }
  return surfacepoint.index;
}

ThinDisk::ThinDisk(Real innerr, Real outerr) :
    innerr(innerr), outerr(outerr) {
}

bool ThinDisk::checkIntersect(const Real &r, const Real &th, const Real &rprev,
    const Real &thprev, SurfacePoint &outSurface) {
  if (r > outerr && rprev > outerr) return false;
  if (r <= innerr && rprev <= innerr) return false;
  if ((th > Pi / 2.0 && thprev < Pi / 2.0)
      || (th < Pi / 2.0 && thprev > Pi / 2.0)) {
    return true;
  }
  return false;
}

int ThinDisk::intersect(const IntegratorData &id,
    const SurfacePoint &surfacepoint, RayHit &hit) {
  hit.r = id.r;
  redshift(hit.r, id.const1, hit.gfactor);
  hit.cosem = id.carter * hit.gfactor / std::sqrt(SQR(hit.r) + epsi3 / hit.r);
  return 512;
}

PlungingRegion::PlungingRegion(Real isco) :
    isco(isco) {
}

bool PlungingRegion::checkIntersect(const Real &r, const Real &th,
    const Real &rprev, const Real &thprev, SurfacePoint &outSurface) {
  if (r > isco && rprev > isco) return false;
  if ((th > Pi / 2.0 && thprev < Pi / 2.0)
      || (th < Pi / 2.0 && thprev > Pi / 2.0)) {
    return true;
  }
  return false;
}
int PlungingRegion::intersect(const IntegratorData &id,
    const SurfacePoint &surfacepoint, RayHit &hit) {
  hit.r = id.r;
  Real met[4][4];
  metric(id.r, id.th, met);
  Real g_tt, g_pp, g_tp;
  g_tt = met[0][0];
  g_pp = met[3][3];
  g_tp = met[0][3];
  Real denom = (g_tt * g_pp - g_tp * g_tp);
  Real ktcalc = -(g_pp + id.b * g_tp) / denom;
  Real kphicalc = (g_tp + id.b * g_tt) / denom;

  Real karray[4] = { ktcalc, id.kr, id.kth, kphicalc };
  redshift_plunge(isco, hit.r, karray, hit.gfactor);
  hit.cosem = id.carter * hit.gfactor / std::sqrt(SQR(hit.r) + epsi3 / hit.r);
  return 600;
}

void Env::addEntity(Entity *entity) {
  this->ents.push_back(entity);
}
bool Env::checkIntersect(const Real &r, const Real &th, const Real &rprev,
    const Real &thprev, SurfacePoint &outsurf, Entity*& hitent) {
  for (Entity *ent : this->ents) {
    if (ent->checkIntersect(r, th, rprev, thprev, outsurf)) {
      hitent = ent;
      return true;
    }
  }
  hitent = nullptr;
  return false;
}
