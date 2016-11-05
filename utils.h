#ifndef __UTILS_H__
#define __UTILS_H__

// UTILS
/** \param v : value to clamp
 * \param a : min edge
 * \param b : max edge
 */
inline float clamp(float v, float a, float b){
	return v<a?a:v>b?b:v;
}

inline float (max)(float a, float b){
  return a>b?a:b;
}

inline float (min)(float a, float b){
  return a<b?a:b;
}

#endif
