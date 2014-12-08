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
ublas::vector<double> calcPixelColor(Ray, ublas::vector<double>, ublas::vector<double>, Material, std::list<Light>, double, double);
#endif /* PIPE_H */
