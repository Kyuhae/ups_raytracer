#ifndef __SCENE_H__
#define __SCENE_H__

#include "vector.h"

// SCENE
typedef struct light_s{
    point3 position;
    color col;
} light;


typedef struct camera_s{
  point3 position;
  vec3 zdir; //! viewdir
  vec3 xdir; //! "left" axis
  vec3 ydir; //! "up axis"
  point3 center;
  float fov; //! field of view in the "xdir" direction
  float aspect; //! width/height of the image
} camera;

// MATERIAL
typedef struct material_s{
    color ks;
    color kd;
    float shininess;  
    float reflect_coef;
} material;


// OBJECT
typedef struct object_s{
  int type;
  material mat;
  //anonymous union of structures to stores object data
  union{
    struct{
      // sphere
      vec3 center;
      float radius;
    };
    struct{
      // plan
      vec3 normal;
      float dist;
    };
    struct{
      // cylindre;
      vec3 base;
      vec3 axis;
      float length;
    };
  };
} object;

camera init_camera(point3 position, point3 at, vec3 up, float fov, float aspect);
object init_cylinder(vec3 base, vec3 axis, float length, float radius, material mat);
object init_sphere(float x, float y, float z, float radius, material mat);
object init_plane(float x, float y, float z, float d, material mat);
#endif
