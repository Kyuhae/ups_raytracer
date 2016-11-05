#ifndef __RAYTRACER_H__
#define __RAYTRACER_H__

#include "ray.h"
#include "scene.h"
#include <stdbool.h>
#include "globals.h"
#include "utils.h"


//! INTERSECTION
//! this structure stores the data corresponding to an intersection, if the intersection is found
typedef struct intersection_s{
  vec3 normal; //! the normal at the intersection point
  point3 position; //! the 3D position of the intersection point
  material mat; //! material of the intersected object
} intersection;

//! render the scene : fill the global variable output_image
void raytrace();


/**
 * Calcul de position d'intersection entre un rayon et un plan.
 * Si l'intersection avec l'objet d'indice passé en paramètre est plus proche que l'intersection actuelle,
 * les informations d'intersection sont mises à jour.
 */
void intersection_plan(ray* current_ray, int ind_obj, intersection* inter);


/**
 * Calcul de position d'intersection entre un rayon et une sphère.
 * Si l'intersection avec l'objet d'indice passé en paramètre est plus proche que l'intersection actuelle,
 *les informations d'intersection sont mises à jour.
 */
void intersection_sphere(ray* current_ray, int ind_obj, intersection* inter);


/**
 * Determine si le point d'intersection donné est eclairé par la lumière donnée.
 * @return: true si ce point est dans l'ombre. false sinon.
 */
bool ombre(intersection* inter, int ind_lum);


/**
 * Calcul de la couleur (directe) du pixel correspondant au rayon passé en paramètre
 * Cette fonction prend en compte le materiau de l'objet, le modèle d'eclairage Blinn-Phong, ainsi que les ombres.
 */
vec3 blinn_phong(ray* current_ray, intersection* inter);


/**
 * Pour un rayon donné, retourne la couleur finale du pixel correspondant
 * Cette fonction ajoute à la couleur directe les reflexions.
 */
vec3 trace(ray* current_ray);


#endif
