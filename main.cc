//Nate Gillard
//CS410 Assignment 3 Fall 2014
//main.cc - entry point for my assignment 3 ray tracer.

#include "objects.h"
#include "pipe.h"
#include "io.h"

int main(int argc, char** argv) {

		
	ublas::vector<float> v (3);
	v(0) = 0;
    v(1) = 3;
	v(2) = -2;

	std::cout<<norm_2(v)<<std::endl;
		

	return 0;
}

