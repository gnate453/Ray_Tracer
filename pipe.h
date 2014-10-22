#ifndef PIPE_H
#define PIPE_H 1

#include <cmath>
#include <algorithm>
#include <sstream>
#include "objects.h"
#include "io.h"

World worldFromString(std::string);
void intersecptRaysandSpheres(World w);

#endif /* PIPE_H */
