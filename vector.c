#include "vector.h"
#include "utils.h"
#include <stdio.h>
#include <math.h>

#ifndef USE_INLINE_MATH
vec3 vector_init(float x, float y, float z){
  vec3 ret;
  ret.x = x;
  ret.y = y;
  ret.z = z;
  return ret;
}

 vec3 vector_minus(vec3 v1, vec3 v2){
  vec3 ret;
  ret.x = v1.x - v2.x;
  ret.y = v1.y - v2.y;
  ret.z = v1.z - v2.z;
  return ret;
}

 vec3 vector_add(vec3 v1, vec3 v2){
  vec3 ret;
  ret.x = v1.x + v2.x;
  ret.y = v1.y + v2.y;
  ret.z = v1.z + v2.z;
  return ret;
}

 vec3 vector_float_mul(float f, vec3 v){
  v.x *= f;
  v.y *= f;
  v.z *= f;
  return v;
}

 vec3 vector_vec_mul(vec3 v1, vec3 v2){
  v1.x *= v2.x;
  v1.y *= v2.y;
  v1.z *= v2.z;
  return v1;
}

 float vector_sqnorm(vec3 v){
  return (v.x*v.x+v.y*v.y+v.z*v.z);
}

 float vector_norm(vec3 v){
    return sqrtf(vector_sqnorm(v));
}

 vec3 vector_normalized(vec3 v){
  float norm = vector_norm(v);
  v.x /= norm;
  v.y /= norm;
  v.z /= norm;
  return v;
}

 vec3 vector_cross(vec3 v1, vec3 v2){
  vec3 ret;
  ret.x = v1.y*v2.z - v1.z*v2.y;
  ret.y = v1.z*v2.x - v1.x*v2.z;
  ret.z = v1.x*v2.y - v1.y*v2.x;
  return ret;
}

 float vector_dot(vec3 v1, vec3 v2){
  return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

 vec3 vector_clamp(vec3 v, float e0, float e1){
  return vector_init(clamp(v.x, e0, e1),
		     clamp(v.y, e0, e1),
		     clamp(v.z, e0, e1));
}

 vec3 vector_reflect(vec3 v, vec3 n){
    return vector_minus(v, vector_float_mul(2.f*vector_dot(n, v), n));
}

 void vector_print(vec3 v){
  printf("v %f %f %f\n", v.x, v.y, v.z);
}

 vec3 vector_min(vec3 v0, vec3 v1){
    return vector_init(min(v0.x, v1.x),
		       min(v0.y, v1.y),
		       min(v0.z, v1.z));
}

vec3 vector_max(vec3 v0, vec3 v1){
    return vector_init(max(v0.x, v1.x),
		       max(v0.y, v1.y),
		       max(v0.z, v1.z));
}
#endif
