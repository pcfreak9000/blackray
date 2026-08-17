#pragma once

#include "def.hpp"
#include <vector>

#define Q1 0
#define Q2 1
#define Q3 2
#define Q4 3

class QuadTree {
public:
  QuadTree(Real x, Real y, Real width, Real height);
  Real check_intersect(Real x1, Real y1, Real x2, Real y2,
      SurfaceElement ** intersect);
  void put_element(SurfaceElement *element);
  void validate();
  size_t size();
private:
  std::vector<QuadTree*> subtrees;
  std::vector<SurfaceElement*> myelements;
  Real x, y, width, height;
  size_t max_elements;
  bool is_leaf;
  size_t level;
//  bool is_root;

  void subdivide();
  bool fits(SurfaceElement *element);
  bool overlaps(Real x0, Real y0, Real x1, Real y1);
  bool fully_inside(Real x0, Real y0, Real x1, Real y1);
};
