#ifndef PIPE_H
#define PIPE_H 1

#include <cmath>
#include <algorithm>
#include <sstream>
#include "objects.h"
#include "io.h"

World worldFromString(std::string);
void castRays(World w);
void intersectRayWithSpheres(Ray, std::list<Sphere>, std::list<Light>, Camera, Image &);
void intersectRayWithPolygons(Ray, std::list<Polygon>, std::list<Light>, Camera, Image &);
#endif /* PIPE_H */
