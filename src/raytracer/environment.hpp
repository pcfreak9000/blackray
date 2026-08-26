#pragma once

#include "def.hpp"
#include "quadtree.hpp"

class Entity {
public:
  virtual int checkIntersect(const Real &r, const Real &th, const Real &rprev,
      const Real &thprev, SurfacePoint &outSurface) = 0;
  virtual int intersect(const IntegratorData &id,
      const SurfacePoint &surfacepoint, RayHit &hit) = 0;

};

class GRMHDDisk : public Entity {
public:
  GRMHDDisk(QuadTree *tree, Real checkr);
  int checkIntersect(const Real &r, const Real &th, const Real &rprev,
      const Real &thprev, SurfacePoint &outSurface) override;
  int intersect(const IntegratorData &id, const SurfacePoint &surfacepoint,
      RayHit &hit) override;
private:
  QuadTree *tree;
  Real checkr;
};

class ThinDisk : public Entity {
public:
  ThinDisk(Real inner, Real outer);
  int checkIntersect(const Real &r, const Real &th, const Real &rprev,
      const Real &thprev, SurfacePoint &outSurface) override;
  int intersect(const IntegratorData &id, const SurfacePoint &surfacepoint,
      RayHit &hit) override;
private:
  Real innerr;
  Real outerr;
};

class PlungingRegion : public Entity {
public:
  PlungingRegion(Real isco);
  int checkIntersect(const Real &r, const Real &th, const Real &rprev,
      const Real &thprev, SurfacePoint &outSurface) override;
  int intersect(const IntegratorData &id, const SurfacePoint &surfacepoint,
      RayHit &hit) override;
private:
  Real isco;
};

class Env {
public:
  Entity* checkIntersect(const Real &r, const Real &th, const Real &rprev,
      const Real &thprev, SurfacePoint &outsurf, int &res);
  void addEntity(Entity *entity);
private:
  std::vector<Entity*> ents;
};

