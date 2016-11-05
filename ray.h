#ifndef __RAY_H__
#define __RAY_H__

#include "vector.h"

// RAY
typedef struct ray_s{
  
  point3 orig; //! start point of the ray
  vec3 dir; //! ray direction, normalized

  float tmax; //! if no intersection computed, should be a large value, else contains the t of the intersection, so do not need to considere intersection if the founded t is more than tmax
  float tmin; //! usefull to store entry point in a aabb
  int depth; //! number of reflection/refraction

} ray;

inline void ray_init(ray *r, vec3 o, vec3 d, float tmin=0, float tmax=1000, int depth=0){
  r->orig = o;
  r->dir = d;
  r->tmin = tmin;
  r->tmax = tmax;
  r->depth = depth;
}

inline vec3 ray_at(ray r, float t){
  return vector_add(vector_float_mul(t, r.dir), r.orig); 
}

#endif
