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

	Image i (2, 1);
	std::cout<<"Hello"<<std::endl;	
	ublas::vector<float> c1 (VECTOR_C);
	c1 (RED) = 134;
	c1 (GREEN) = 0.0;
	c1 (BLUE) = 0.0;
	c1 (ALPHA) = 1.0;
	std::cout<<"Hello"<<std::endl;	
	i.setPixelColor(0, 0, c1);
	i.setPixelColor(1, 0, c1);
	
	std::cout<<"Hello"<<std::endl;	
	i.setPixelDepth(1, 0, 345.7);
	i.setPixelDepth(0, 0, 1.0);

	std::cout<<i.getPixelColor(0,0)<<std::endl;
	std::cout<<i.getPixelColor(1,0)<<std::endl;
	std::cout<<i.getPixelDepth(0,0)<<std::endl;
	std::cout<<i.getPixelDepth(1,0)<<std::endl;
	return 0;
}

