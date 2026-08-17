#include <iostream>

#include "quadtree.hpp"

Real checkIntersect(Real x1, Real y1, Real x2,
    Real y2, Real x3, Real y3, Real x4,
    Real y4) {
  Real axl = std::min(x1, x2);
  Real axh = std::max(x1, x2);
  Real ayl = std::min(y1, y2);
  Real ayh = std::max(y1, y2);
  Real bxl = std::min(x3, x4);
  Real bxh = std::max(x3, x4);
  Real byl = std::min(y3, y4);
  Real byh = std::max(y3, y4);
  if (!(axl < bxh && axh > bxl && ayl < byh && ayh > byl)) {
    return NO_INTERSECT;
  }
  Real num = (x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4);
  Real denum = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
  Real t = num / denum;
  if (t >= 0 && t <= 1) {
    return t;
  }
  return NO_INTERSECT;
}

QuadTree::QuadTree(Real x, Real y, Real width, Real height) :
    x(x), y(y), width(width), height(height), max_elements(4), is_leaf(true),level(0) {
}

Real QuadTree::check_intersect(Real x1, Real y1, Real x2, Real y2, SurfaceElement** out){

  if(!is_leaf) {
    for(int i=0; i<4; i++){
      if(subtrees[i]->overlaps(x1,y1,x2,y2)){
        Real ressub = subtrees[i]->check_intersect(x1,y1,x2,y2,out);
        if(ressub != NO_INTERSECT){
          return ressub;
        }
        //this needs some benchmarking to check if this actually improves things:
        if(subtrees[i]->fully_inside(x1,y1,x2,y2)){
          break;
        }
      }
    }
  }
  for(int i=0; i<myelements.size(); i++) {
        SurfaceElement* elem = myelements[i];
        if (elem->sp0->y==0.0 && elem->sp1->y==0.0){
          //if this raytracer is generalized, this does not really belong here. But we ignore segments with y:=0 of the accretion disk.
          continue;
        }
        Real res = checkIntersect(elem->sp0->x,elem->sp0->y,elem->sp1->x,elem->sp1->y,x1,y1,x2,y2);
        if(res != NO_INTERSECT) {
          *out = elem;
          return res;
        }
    }

  return NO_INTERSECT;
}

size_t QuadTree::size(){
  if(is_leaf) return myelements.size();
  size_t s = myelements.size();
  for(int i=0; i<4; i++){
    s += subtrees[i]->size();
  }
  return s;
}

void QuadTree::validate(){
  //std::cout << "Validating..." << std::endl;
  for(SurfaceElement* se : myelements){
    if(!se) std::cout << level << " " << myelements.size() << "--- uh oh ---: " << se << std::endl;
    if(!se->sp0) std::cout << level << " " << myelements.size() << " oh no 0: " << se->sp0 << std::endl;
    if(!se->sp1) std::cout << level << " " << myelements.size() << " oh no 1: " << se->sp1 << std::endl;
  }
  if(!is_leaf){
    for(QuadTree* qt : subtrees){
      qt->validate();
    }
  }
}

bool QuadTree::fits(SurfaceElement *element) {
  Real maxx = std::max(element->sp0->x, element->sp1->x);
  Real minx = std::min(element->sp0->x, element->sp1->x);
  Real maxy = std::max(element->sp0->y, element->sp1->y);
  Real miny = std::min(element->sp0->y, element->sp1->y);

  return maxx < x+width && minx >= x && maxy < y+height && miny >= y;
}

bool QuadTree::overlaps(Real x0, Real y0, Real x1, Real y1) {
  Real maxx = std::max(x0, x1);
  Real minx = std::min(x0, x1);
  Real maxy = std::max(y0, y1);
  Real miny = std::min(y0, y1);
  return maxx >= x && minx < x+width && maxy >= y && miny < y+height;
}

bool QuadTree::fully_inside(Real x0, Real y0, Real x1, Real y1) {
  Real maxx = std::max(x0, x1);
  Real minx = std::min(x0, x1);
  Real maxy = std::max(y0, y1);
  Real miny = std::min(y0, y1);

  return maxx < x+width && minx >= x && maxy < y+height && miny >= y;
}
void QuadTree::subdivide() {
  Real w_2 = width/2.0;
  Real h_2 = height/2.0;
  QuadTree* q1 = new QuadTree(x+w_2, y+h_2, w_2, h_2);
  QuadTree* q2 = new QuadTree(x,     y+h_2, w_2, h_2);
  QuadTree* q3 = new QuadTree(x,     y,     w_2, h_2);
  QuadTree* q4 = new QuadTree(x+w_2, y,     w_2, h_2);
  q1->level = this->level+1;
  q2->level = this->level+1;
  q3->level = this->level+1;
  q4->level = this->level+1;
  is_leaf = false;
  subtrees.push_back(q1);
  subtrees.push_back(q2);
  subtrees.push_back(q3);
  subtrees.push_back(q4);
  for(auto vecit = myelements.begin(); vecit != myelements.end(); vecit++){
    SurfaceElement* se = *vecit;
    for(int i=0; i<4; i++){
      if(subtrees[i]->fits(se)){
        subtrees[i]->put_element(se);
        myelements.erase(vecit);
        vecit--;
        break;
      }
    }
  }
}

void QuadTree::put_element(SurfaceElement *element){
  if(is_leaf) {
    myelements.push_back(element);
    if(myelements.size() > max_elements){
      subdivide();
    }
  } else {
    for(int i=0; i<4; i++){
      if(subtrees[i]->fits(element)){
        subtrees[i]->put_element(element);
        break;
      }
      if(i==3){
        myelements.push_back(element);
      }
    }
  }
}
