#include "raytracer.h"
#include <math.h>


void intersection_plan(ray* current_ray, int ind_obj, intersection* inter){ 
  float t, tmp;

  tmp = vector_dot(scene[ind_obj].normal, current_ray->orig) + scene[ind_obj].dist;
  t = -tmp / vector_dot(current_ray->dir, scene[ind_obj].normal);
  
  if (t > 0 && t < current_ray->tmax) {
    current_ray->tmax = t;
    inter->normal = scene[ind_obj].normal;
    inter->position = ray_at(*current_ray, t);
    inter->mat = scene[ind_obj].mat;
  }
}



void intersection_sphere(ray* current_ray, int ind_obj, intersection* inter){
  float t1, t2, a, b, c, d;
  vec3 vecCO;

  //Resolution equation du 2nd degré pour avoir t
  vecCO = vector_minus(current_ray->orig, scene[ind_obj].center);
  a = vector_sqnorm(current_ray->dir);
  b = 2 * vector_dot(current_ray->dir, vecCO);
  c = vector_sqnorm(vecCO) - scene[ind_obj].radius * scene[ind_obj].radius  ;
  
  d = (b*b) - 4 * a * c;
  
  //traitement différents cas. Sortie si aucune solution convenable.
  if (d < 0) return;
  if (d == 0) {
    t1 = -b / (2 * a);
    if (t1 > 0 && t1 < current_ray->tmax) {
      current_ray->tmax = t1;
      inter->position = ray_at(*current_ray, t1);
      inter->normal = vector_normalized(vector_minus(inter->position, scene[ind_obj].center));
      inter->mat = scene[ind_obj].mat;
    }
  }
  else
    if (d > 0) {
      t1 = (- b - sqrt(d)) / (2 * a);
      t2 = (- b + sqrt(d)) / (2 * a);
      if (t1 > 0 && t1 < current_ray->tmax && t2 > 0 && t2 < current_ray->tmax) {
        current_ray->tmax = fmin(t1, t2);
        inter->position = ray_at(*current_ray, fmin(t1, t2));
        inter->normal = vector_normalized(vector_minus(inter->position, scene[ind_obj].center));
        inter->mat = scene[ind_obj].mat;
      }
      if (t1 > 0 && t1 < current_ray->tmax && (t2 < 0 || t2 >current_ray->tmax)) {
        current_ray->tmax = t1;
        inter->position = ray_at(*current_ray, t1);
        inter->normal = vector_normalized(vector_minus(inter->position, scene[ind_obj].center));
        inter->mat = scene[ind_obj].mat;
      }
      if (t2 > 0 && t2 < current_ray->tmax && (t1 < 0 || t1 >current_ray->tmax)) {
        current_ray->tmax = t2;
        inter->position = ray_at(*current_ray, t2);
        inter->normal = vector_normalized(vector_minus(inter->position, scene[ind_obj].center));
        inter->mat = scene[ind_obj].mat;
      }
    }
}


bool ombre(intersection* inter, int ind_lum) {
  intersection interOmbre;
  ray current_ray;
  bool inter_found;
  vec3 dir, pos;
  
  //Creation du vecteur Intersection ---> Lumière
  dir = vector_minus(lights[ind_lum].position, inter->position);
  pos = vector_add(inter->position, vector_float_mul(EPSILON, vector_normalized(dir)));
  
  ray_init(&current_ray, pos, dir );
  current_ray.tmax = 1;
  
  //Determiner l'intersection la plus proche
  inter_found = false;
  for (int k = 0; k < object_count && !inter_found; k++) {
    if (scene[k].type == PLANE){
      intersection_plan(&current_ray, k, &interOmbre);
    }
    else if (scene[k].type == SPHERE)
      intersection_sphere(&current_ray, k, &interOmbre);
    if (!inter_found && current_ray.tmax != 1) 
      inter_found = true;
   }
   
  return inter_found; 
}



vec3 blinn_phong(ray* current_ray, intersection* inter) {
  vec3 color = vector_init(0, 0, 0);
  vec3 nul = vector_init(0, 0, 0);
  vec3 dir_light, dir_view, vec_h, i_TIMES_nDOTl, kd_OVER_pi, ks_TIMES_coef;
  float hDOTn_POWs;
  
  for (int i = 0; i < num_lights; i++) {
    if (!ombre(inter, i)) {
      dir_light = vector_normalized(vector_minus(lights[i].position, inter->position));
      dir_view = vector_normalized(vector_minus(nul, current_ray->dir));
      vec_h = vector_normalized(vector_add(dir_light, dir_view));
      
      i_TIMES_nDOTl = vector_float_mul(clamp(vector_dot(inter->normal, dir_light),0.0, 1.0), lights[i].col);
      kd_OVER_pi = vector_float_mul((1 / M_PI), inter->mat.kd);
      hDOTn_POWs = pow(clamp(vector_dot(vec_h, inter->normal), 0.0, 1.0), inter->mat.shininess);
      ks_TIMES_coef = vector_float_mul(( (inter->mat.shininess + 8) / (8 * M_PI)), inter->mat.ks);
      
      color = vector_add(color, vector_vec_mul(i_TIMES_nDOTl, vector_add(kd_OVER_pi, vector_float_mul(hDOTn_POWs, ks_TIMES_coef))));
    }
  }
  return color;   
}

vec3 trace(ray* current_ray) {
  intersection inter;
  bool inter_found;
  vec3 dir, pos, color;
  vec3 sky_blue = vector_init(0.2, 0.2, 0.6);
  ray reflected_ray;
  
  //Determiner l'intersection la plus proche
  inter_found = false;
  for (int k = 0; k < object_count; k++) {
    if (scene[k].type == PLANE)
      intersection_plan(current_ray, k, &inter);
    else if (scene[k].type == SPHERE)
      intersection_sphere(current_ray, k, &inter);
    if (!inter_found && current_ray->tmax != 1000) 
      inter_found = true;
  }

  //Initialisation du rayon refléchi
  dir = vector_reflect(current_ray->dir, inter.normal);
  pos = vector_add(inter.position, vector_float_mul(EPSILON, vector_normalized(dir)));
  
  ray_init(&reflected_ray, pos, dir);
  reflected_ray.depth = current_ray->depth + 1;
  
  //Definir la couleur de ce pixel, avec gestion des reflets
  if (inter_found) {
    if (current_ray->depth < 5)
      return (vector_add(blinn_phong(current_ray, &inter), vector_vec_mul(inter.mat.ks, trace(&reflected_ray)))); 
    else
      return vector_init(0, 0, 0);
  }
  else{    
    return sky_blue;
  } 
}


void raytrace() {
  ray current_ray;
  vec3 xr, yr, zr, direction, color;
  vec3 sky_blue = vector_init(0.2, 0.2, 0.6);
  float xcoef, ycoef;
  bool inter_found;
  intersection inter;
  
  for (int j = 0; j < HEIGHT; j++) {
    for (int i = 0; i < WIDTH; i++) {
      
      //determiner l'equation du rayon
      xr = theCamera.xdir;
      yr = vector_float_mul(1 / theCamera.aspect, theCamera.ydir);
      zr = vector_float_mul(1 / tan(0.5 * ((theCamera.fov /180)* M_PI)), theCamera.zdir);
      
      xcoef = (i + 0.5 - (float)(WIDTH / 2)) / ((float)(WIDTH / 2));
      ycoef = (j + 0.5 - (float)(HEIGHT / 2)) / (float)((HEIGHT / 2));

      //Initialisation du rayon
      direction = vector_add(zr, vector_add(vector_float_mul(xcoef, xr), vector_float_mul(ycoef, yr)));
      ray_init(&current_ray, theCamera.position, direction);
      
      //Calcul de la couleur du pixel associé à ce rayon
      color = trace(&current_ray);
      
      output_image[j * WIDTH + i] = color;
      
      if (i == WIDTH - 1)
        printf("[%d%] \r", (i * j * 100) / (WIDTH * HEIGHT));
       
    }
  }
}
