#ifndef PIPE_H
#define PIPE_H 1

#include <cmath>
#include <algorithm>
#include <sstream>
#include "objects.h"
#include "io.h"

float smoothDecimal(float, int);
ublas::vector<float> smoothColor(ublas::vector<float>);
World worldFromString(std::string);
void castRays(World w);
ublas::vector<float> rayTrace(Ray, World*, int, bool);
Intersection* intersectRayWithSpheres(Ray, std::list<Sphere>);
Intersection* intersectRayWithPolygons(Ray, std::list<Polygon>);
ublas::vector<float> calcSpecularColor(Material, ublas::vector<float>, float);
ublas::vector<float> calcDiffuseColor(Material, ublas::vector<float>, float);
ublas::vector<float> calcAmbientColor(Material);
#endif /* PIPE_H */
