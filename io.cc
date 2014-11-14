#include "io.h"

std::string readInputFiles(char** argv) {

	std::string tmp, mtlFileName, objFileName, materials, vertices, objects, commands, world;

	//open the object file for reading.
	objFileName = argv[1];
	std::ifstream objFile(objFileName.c_str());
	
	//get the material file name.
	std::getline(objFile, tmp);
	tmp = tmp.substr(LENGTH_7 , (tmp.size() - 1));
	
	if (objFileName.find("\\/") != std::string::npos) {
		mtlFileName = objFileName.substr(BEGIN, objFileName.find_last_of("\\/"));
		mtlFileName = mtlFileName + tmp;
	}
	else {
		mtlFileName = tmp;
	}
	
	//get the rest of the object file.
	while (std::getline(objFile, tmp))
	{
		if (*(tmp.begin()) == 'v')
			vertices += tmp + "\n";
		else if (*(tmp.begin()) == 'f')
			objects += tmp + "\n";
		else if (*(tmp.begin()) == 'g')
			objects += tmp + "\n";
		else if (*(tmp.begin()) == 's')
			objects += tmp + "\n";
		else if (*(tmp.begin()) == 'u')
			objects += tmp + "\n";
	}
	objFile.close();	

	//open the material file for reading.
	std::ifstream mtlFile(mtlFileName.c_str());
	//get lines from materials file.
	while (std::getline(mtlFile, tmp))
	{
		std::string check = tmp.substr(BEGIN, LENGTH_2);
		if ( check.compare("n1") != 0 &&
			 check.compare("Tr") != 0 &&
			 check.compare("Kr") != 0 )		
			materials += tmp + "\n";
	}
	mtlFile.close();

	//open the command file for reading.
	std::ifstream commFile(argv[2]);

	//get lines from command file.
	while (std::getline(commFile, tmp))
	{
		commands += tmp + "\n";
	}
	commFile.close();
	
	world = vertices + materials + objects + commands;

	//std::cout<<world<<std::endl;
	//std::cout<<"parse 1 done"<<std::endl;

	return world;
}

void outputImages(std::list<Image> imgs) {
	
	for (std::list<Image>::iterator iit = imgs.begin();
			iit != imgs.end(); ++iit)
	{
		std::string nameC = iit->getName() + "_color.ppm";
		std::string nameD = iit->getName() + "_depth.ppm";
		int w = iit->getWidth();
		int h = iit->getHeight();
		
		std::ofstream colorOut(nameC.c_str(), std::ofstream::out );
		std::ofstream depthOut(nameD.c_str(), std::ofstream::out );
		
		//std::cout<<w<<" "<<h<<std::endl;

		colorOut<<"P3 "<<w<<" "<<h<<" "<<(COLOR_MAX+1)<<std::endl;
		depthOut<<"P3 "<<w<<" "<<h<<" "<<(COLOR_MAX+1)<<std::endl;

		for (int y = 0; y < h; ++y) 
		{
			for (int x = 0; x < w; ++x)
			{
			 
				colorOut<<iit->getPixelRed(x, y)<<" "
						<<iit->getPixelGreen(x, y)<<" "
						<<iit->getPixelBlue(x, y)<<" ";

				depthOut<<(int) floor(iit->getPixelDepth(x, y))<<" "
						<<(int) floor(iit->getPixelDepth(x, y))<<" "
						<<(int) floor(iit->getPixelDepth(x, y))<<" ";
			}
			colorOut<<"\n";
			depthOut<<"\n";
		}

		colorOut.close();
		depthOut.close();
	}	

}
