//Nate Gillard
//CS410 Assignment 3 Fall 2014
//main.cc - entry point for my assignment 3 ray tracer.

#include "objects.h"
#include "pipe.h"
#include "io.h"

int main(int argc, char** argv) {		
	
	if (argc < 3 || argc < 3)
	{
		std::cerr<<"Usage: \n HW5 <.obj> <command>"<<std::endl;
		return 1;
	}

	std::string data = readInputFiles(argv);
	
	World w = worldFromString(data);

	castRays(w);

	return 0;
}

