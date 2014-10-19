//Nate Gillard
//CS410 Assignment 3 Fall 2014
//main.cc - entry point for my assignment 3 ray tracer.

#include "objects.h"
#include "pipe.h"
#include "io.h"

int main(int argc, char** argv) {		

	ublas::vector<float> prp (3);
	ublas::vector<float> pixel (3);

	prp (0) = 0;
	prp (1) = 0;
	prp (2) = 0;
	pixel (0) = 0;
	pixel (1) = 3;
	pixel (2) = -2;

	Ray r(prp, pixel); 

	std::cout<<r.getPixel()<<std::endl;
	std::cout<<r.getPRP()<<std::endl;
	std::cout<<r.unitVector()<<std::endl;
	
	return 0;
}

