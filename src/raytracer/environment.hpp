#pragma once

class Entity;
class Env;

#include "def.hpp"
#include "quadtree.hpp"
#include "raytracingnew.hpp"

class Entity {
public:
  virtual int checkIntersect(const Real &r, const Real &th, const Real &rprev,
      const Real &thprev) = 0;
  virtual Real getMaxRadius() = 0;

  virtual int calculateRedshift(const InitialCondition* ic, const RayHit &hit, Real& gfactor, Real& cosem);
};

class GRMHDDisk : public Entity {
public:
  GRMHDDisk(QuadTree *tree, Real checkr);
  int checkIntersect(const Real &r, const Real &th, const Real &rprev,
      const Real &thprev) override;
  Real getMaxRadius() override;
  int calculateRedshift(const InitialCondition* ic, const RayHit &hit, Real& gfactor, Real& cosem) override;
private:
  QuadTree *tree;
  Real checkr;
};

class ThinDisk : public Entity {
public:
  ThinDisk(Real inner, Real outer);
  int checkIntersect(const Real &r, const Real &th, const Real &rprev,
      const Real &thprev) override;
  Real getMaxRadius() override;
  int calculateRedshift(const InitialCondition* ic, const RayHit &hit, Real& gfactor, Real& cosem) override;
private:
  Real innerr;
  Real outerr;
};

class PlungingRegion : public Entity {
public:
  PlungingRegion(Real isco);
  int checkIntersect(const Real &r, const Real &th, const Real &rprev,
      const Real &thprev) override;
  Real getMaxRadius() override;
  int calculateRedshift(const InitialCondition* ic, const RayHit &hit, Real& gfactor, Real& cosem) override;
private:
  Real isco;
};

class Env {
public:
  int checkIntersect(const Real &r, const Real &th, const Real &rprev,
      const Real &thprev, Entity*& hitent);
  void addEntity(std::unique_ptr<Entity> ptr);
private:
  std::vector<std::unique_ptr<Entity>> ents;
  Real maxr = 0;
};

