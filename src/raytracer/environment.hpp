#pragma once

#include "def.hpp"
#include "quadtree.hpp"

class Entity {
public:
  virtual bool checkIntersect(const Real &r, const Real &th, const Real &rprev,
      const Real &thprev, SurfacePoint &outSurface) = 0;
  virtual int intersect(const IntegratorData &id,
      const SurfacePoint &surfacepoint, RayHit &hit) = 0;
  virtual Real getMaxRadius() = 0;
};

class GRMHDDisk : public Entity {
public:
  GRMHDDisk(QuadTree *tree, Real checkr);
  bool checkIntersect(const Real &r, const Real &th, const Real &rprev,
      const Real &thprev, SurfacePoint &outSurface) override;
  int intersect(const IntegratorData &id, const SurfacePoint &surfacepoint,
      RayHit &hit) override;
  Real getMaxRadius() override;
private:
  QuadTree *tree;
  Real checkr;
};

class ThinDisk : public Entity {
public:
  ThinDisk(Real inner, Real outer);
  bool checkIntersect(const Real &r, const Real &th, const Real &rprev,
      const Real &thprev, SurfacePoint &outSurface) override;
  int intersect(const IntegratorData &id, const SurfacePoint &surfacepoint,
      RayHit &hit) override;
  Real getMaxRadius() override;
private:
  Real innerr;
  Real outerr;
};

class PlungingRegion : public Entity {
public:
  PlungingRegion(Real isco);
  bool checkIntersect(const Real &r, const Real &th, const Real &rprev,
      const Real &thprev, SurfacePoint &outSurface) override;
  int intersect(const IntegratorData &id, const SurfacePoint &surfacepoint,
      RayHit &hit) override;
  Real getMaxRadius() override;
private:
  Real isco;
};

class Env {
public:
  bool checkIntersect(const Real &r, const Real &th, const Real &rprev,
      const Real &thprev, SurfacePoint &outsurf, Entity*& hitent);
  void addEntity(Entity *entity);
private:
  std::vector<Entity*> ents;
  Real maxr = 0;
};

