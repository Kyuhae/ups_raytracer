#include "scene.h"
#include "globals.h"

object init_cylinder(vec3 base, vec3 axis, float length, float radius, material mat){
  object ret;
  ret.type = CYLINDER;
  ret.base = base;
  ret.axis = vector_normalized(axis);
  ret.radius = radius;
  ret.length = length;  
  ret.mat = mat;
  return ret;
}

object init_sphere(float x, float y, float z, float radius, material mat){
  object ret;
  ret.type = SPHERE;
  ret.center = vector_init(x, y, z);
  ret.radius = radius;
  ret.mat = mat;
  return ret;
}

object init_plane(float x, float y, float z, float d, material mat){
  object ret;
  ret.type = PLANE;
  ret.normal = vector_normalized(vector_init(x, y, z));
  ret.dist = d;
  ret.mat = mat;
  return ret;
}



camera init_camera(point3 position, point3 at, vec3 up, float fov, float aspect){
    camera cam;
    cam.fov = fov;
    cam.aspect = aspect;
    cam.position = position;

    cam.zdir = vector_normalized(vector_minus(at, position));
    cam.xdir = vector_normalized(vector_cross(up, cam.zdir));
    cam.ydir = vector_normalized(vector_cross(cam.zdir, cam.xdir));
    cam.center = vector_float_mul( 1.f / tanf ((cam.fov * M_PI / 180.f) * 0.5f) , cam.zdir);
    return cam;
}
