#ifndef IO_H
#define IO_H 1

#include <string>
#include <fstream>

#include "objects.h"

std::string readInputFiles(char**);
void outputImages(std::list<Image>);

#endif	/* IO_H */
