#ifndef __VECTOR_H__
#define __VECTOR_H__


#include "utils.h"
#include <stdio.h>
#include <math.h>


// VECTORS
typedef struct vector_s{

  // first coordinate access as .value[0] .x or .r
    union{
	float value[3];
	struct{
            float x;
            float y;
            float z;
        };
    };
} vec3;


typedef vec3 color;
typedef vec3 point3;

#ifndef USE_INLINE_MATH
vec3 vector_init(float x, float y, float z);
vec3 vector_minus(vec3 v1, vec3 v2);
vec3 vector_add(vec3 v1, vec3 v2);
vec3 vector_float_mul(float f, vec3 v);
vec3 vector_vec_mul(vec3 v1, vec3 v2);
float vector_sqnorm(vec3 v);
float vector_norm(vec3 v);
vec3 vector_normalized(vec3 v);
vec3 vector_cross(vec3 v1, vec3 v2);
float vector_dot(vec3 v1, vec3 v2);
vec3 vector_clamp(vec3 v, float e0, float e1);
vec3 vector_reflect(vec3 v, vec3 n);
void vector_print(vec3 v);
 //! \brief min component wise
vec3 vector_min(vec3 v0, vec3 v1);
//! \brief max component wise
vec3 vector_max(vec3 v0, vec3 v1);
#else
inline vec3 vector_init(const float x, const float y, const float z){
  vec3 ret;
  ret.x = x;   ret.y = y; ret.z = z;
  return ret;
}

inline vec3 vector_minus(const vec3 v1, const vec3 v2){
  return vector_init( v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

inline vec3 vector_add(const vec3 v1, const vec3 v2){
  return vector_init(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

inline vec3 vector_float_mul(const float f, const vec3 v){
  return vector_init( v.x * f,   v.y * f,  v.z * f);
}

inline vec3 vector_vec_mul(const vec3 v1, const vec3 v2){
  return vector_init( v1.x * v2.x,  v1.y * v2.y,  v1.z * v2.z);
}

inline float vector_sqnorm(const vec3 v){
  return (v.x*v.x+v.y*v.y+v.z*v.z);
}

inline float vector_norm(const vec3 v){
    return sqrtf(vector_sqnorm(v));
}

inline vec3 vector_normalized(const vec3 v){
  float norm = vector_norm(v);
  return vector_init(v.x / norm,   v.y / norm,  v.z / norm);
}

//! \todo passage par pointeur
inline vec3 vector_cross(const vec3 v1, const vec3 v2){
  return vector_init(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x);
}

inline float vector_dot(const vec3 v1, const vec3 v2){
  return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

inline vec3 vector_clamp(const vec3 v, const float e0, const float e1){
    return vector_init(clamp(v.x, e0, e1), clamp(v.y, e0, e1), clamp(v.z, e0, e1));
}

inline vec3 vector_reflect(const vec3 v, const vec3 n){
    return vector_minus(v, vector_float_mul(2.f*vector_dot(n, v), n));
}

inline void vector_print(const vec3 v){
    printf("v %f %f %f\n", v.x, v.y, v.z);
}

inline vec3 vector_min(const vec3 v0, const vec3 v1){
    return vector_init(min(v0.x, v1.x), min(v0.y, v1.y), min(v0.z, v1.z));
}

inline vec3 vector_max(vec3 v0, vec3 v1){
    return vector_init(max(v0.x, v1.x), max(v0.y, v1.y), max(v0.z, v1.z));
}
#endif
#endif
