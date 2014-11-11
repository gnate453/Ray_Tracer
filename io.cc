#include "io.h"

std::string readInputFiles(char** argv) {

	std::string tmp, mtlFileName, objFileName, materials, vertices, objects, commands, world;

	//open the object file for reading.
	objFileName = argv[1];
	std::ifstream objFile(objFileName);
	
	//get the material file name.
	std::getline(objFile, tmp);
	mtlFileName = objFileName.substr(0, objFileName.find_last_of("\\/");
	mtlFileName = mtlFileName + tmp;

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
	std::ifstream mtlFile(mtlFileName);

	//get lines from materials file.
	while (std::getline(mtlFile, tmp))
	{
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
