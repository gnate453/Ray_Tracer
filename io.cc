#include "io.h"

std::string readInputFiles(char** argv) {

	std::ifstream objIn(argv[1]);
	std::ifstream cmdIn(argv[2]);

	std::string objSTR;
	std::string tmp;

	while (std::getline(objIn, tmp))
	{
		objSTR += tmp + "\n";
	}

	while (std::getline(cmdIn, tmp))
	{
		objSTR += tmp + "\n";
	}
	
	return objSTR;
}
