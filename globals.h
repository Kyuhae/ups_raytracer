#ifndef __GLOBALS_H__
#define __GLOBALS_H__

#include "scene.h"
#include "vector.h"
#include "raytracer.h"

#ifndef M_PI 
#define M_PI 3.14159
#endif

#define EPSILON 0.0001f
#define WIDTH (1280)
#define HEIGHT (720)

#define CYLINDER 3
#define PLANE 2
#define SPHERE 1

#define SCENE_SIZE (1024)

extern light lights[3];
extern int num_lights;

extern object scene[SCENE_SIZE];
extern int object_count;
extern color output_image[WIDTH*HEIGHT];
extern camera theCamera;


#endif
