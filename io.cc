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
		
		std::cout<<w<<" "<<h<<std::endl;

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
