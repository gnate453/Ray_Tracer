#ifndef PIPE_H
#define PIPE_H 1

#include <cmath>
#include <algorithm>
#include <sstream>
#include "objects.h"
#include "io.h"

World worldFromString(std::string);
void castRays(World w);
Intersection intersectRayWithSpheres(Ray, std::list<Sphere>, Camera);
Intersection intersectRayWithPolygons(Ray, std::list<Polygon>, Camera);
ublas::vector<float> calcPixelColor(Ray, ublas::vector<float>, ublas::vector<float>, Material, std::list<Light>, float, float);
#endif /* PIPE_H */
