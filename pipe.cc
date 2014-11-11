#include "pipe.h"

World worldFromString(std::string data) {
	
	World newWorld;
	std::stringstream d(data);
	std::string line, entry, name;
	
	while (std::getline(d, line))
	{
		std::stringstream l(line);
		std::getline(l, entry, ' ');

		if (entry.compare("s") == 0)
		{
			std::getline(l, name, ' ');
			ublas::vector<float> w (VECTOR_3D);
			l>>w (X);
			l>>w (Y);
			l>>w (Z);
			float r;
			l>>r;
			ublas::vector<float> c (VECTOR_C);
			l>>c (RED);
			l>>c (GREEN);
			l>>c (BLUE);
			c (ALPHA) = 0;
			//std::cout<<name<<" "<<w<<" "<<c<<" "<<r<<std::endl;
			Sphere sphere(name, w, c, r);
			newWorld.addSphere(sphere);
		}
		else if (entry.compare("c") == 0)
		{
			std::getline(l, name, ' ');
			ublas::vector<float> prp (VECTOR_3D);
			l>>prp (X);
			l>>prp (Y);
			l>>prp (Z);
			ublas::vector<float> vpn (VECTOR_3D);
			l>>vpn (X);
			l>>vpn (Y);
			l>>vpn (Z);
			float near, far;
			l>>near;
			l>>far;
			//put the vpn at z= -near;
			vpn (Z) = -near; 
			Camera camera(name, prp, vpn, near, far);
			newWorld.addCamera(camera);
	
		}
		else if (entry.compare("r") == 0)
		{
			getline(l, name, ' ');
			int width, height, depth;
			l>>width;
			l>>height;
			l>>depth;
			Scene scene(name, width, height, depth);
			newWorld.addScene(scene);
		}
	}

	return newWorld;
}

void castRays(World w) {

	std::list<Sphere> spheres = w.getSpheres();
	std::list<Camera> cameras = w.getCameras();
	std::list<Scene> scenes = w.getScenes();
	std::list<Image> imgs;
	float rsqd, csqd;

	
	for (std::list<Scene>::iterator sit = scenes.begin();
			sit != scenes.end(); ++sit)
	{
		for (std::list<Camera>::iterator cit = cameras.begin();
				cit != cameras.end(); ++cit)
		{
			//get r^2 and c^2 for each sphere.
			for (std::list<Sphere>::iterator oit = spheres.begin();
					oit != spheres.end(); ++oit)
			{
				//radius squared
				rsqd = (oit->getRadius()*oit->getRadius());
				oit->setRadiusSquared(rsqd);
			
				//csqd = norm_2(oit->getOrigin() -  cit->getPRP()) * norm_2(oit->getOrigin() -  cit->getPRP());
				//oit->setDistanceToPRPSquared(csqd);

				//std::cout<<oit->getName()<<" rsqd: "<<rsqd<<" "<<"csqd: "<<csqd<<std::endl;	
			}//end for each sphere
	
			Image img(sit->getName()+"_"+cit->getName(), sit->getWidth(), sit->getHeight());
			for (int i = 0; i < sit->getHeight(); ++i)
			{
				float y = (2.0/(sit->getHeight()-1))*i - 1 ;
				for (int j = 0; j < sit->getWidth(); ++j)
				{
					Ray tmpR(cit->getPRP()) ;
					float x = (2.0/(sit->getWidth()-1))*j - 1;
					tmpR.setX(j);
					tmpR.setY(i);
					//std::cout<<x<<","<<y<<" ";
					//std::cout<<j<<","<<i<<" ";
					ublas::vector<float> rv (VECTOR_3D);
					rv (X) = x;
					rv (Y) = y;
					rv (Z) = -cit->getNearClip();
					tmpR.setPixel(rv);
					intersectRayWithObjects(tmpR, spheres, *cit, img);
				}
				//std::cout<<std::endl;
			}	
					
			imgs.push_back(img);

		}//end for each camera
	}//end for each scene


	outputImages(imgs);
}

void intersectRayWithObjects(Ray ray, std::list<Sphere> objs, Camera c, Image &img) {

	ublas::vector<float> U = ray.unitVector();
	//std::cout<<ray.getX()<<","<<ray.getY()<<" U: "<<U<<std::endl;
	int pr = 0;
	int pg = 0;
	int pb = 0;
	for (std::list<Sphere>::iterator oit = objs.begin(); oit != objs.end(); ++oit)
	{
		float v, csqd, dsqd;
		//std::cout<<"Object: "<<oit->getName()<<std::endl;
		//c squared (sphere origin -camera prp)
		if (ray.getPixel() (Z) > oit->getOrigin() (Z)) {
		//	std::cout<<"sphere past near clip "<<std::endl;
			csqd = norm_2(oit->getOrigin() - ray.getPixel()) *  norm_2(oit->getOrigin() - ray.getPixel());
			v = inner_prod(oit->getOrigin() - ray.getPixel(), U);
			dsqd = oit->getRadiusSquared() - ( csqd - (v * v) );
		}
		else {
		//	std::cout<<"sphere before near clip "<<std::endl;
			csqd = norm_2(ray.getPixel() - oit->getOrigin()) *  norm_2(ray.getPixel() - oit->getOrigin());
			v = inner_prod(ray.getPixel() - oit->getOrigin() , U);
			dsqd = ( csqd - (v * v) ) - oit->getRadiusSquared();
		}
		oit->setDistanceToPRPSquared(csqd);
		//std::cout<<" r^2: "<<oit->getRadiusSquared()<<" c^2: "<<oit->getDistanceToPRPSquared()<<" v^2: "<<v*v<<std::endl;
		if (dsqd >= 0)
		{
			float d = std::sqrt(dsqd);
		//	std::cout<<"  d^2: "<<dsqd<<" d: "<<d<<" v: "<<v<<"\n"<<std::endl;
			//normal to sphere
			ublas::vector<float> S = ray.paraPos(v-d);
			if ( S(Z) < -c.getNearClip()  && S(Z) > -c.getFarClip()) { 
				//std::cout<<"Object: "<<oit->getName()<<" S(Z): "<<S(Z)<<std::endl;
				//std::cout<<"near: "<<-cit->getNearClip()<<" far: "<<-cit->getFarClip()<<std::endl;
				ublas::vector<float> N = oit->getOrigin() - S;
				N = (1/norm_2(N)) * N; 
						
				float k = inner_prod(ray.unitVector() , N);
				//std::cout<<"  Q: "<<Q<<" k: "<<k<<std::endl;
						
				int r = (int) floor(oit->getColor() (RED) * std::abs(k) + AMB_LIGHT);
				if (r < COLOR_MAX)
					pr = r;
				else 
					pr = COLOR_MAX;

				int g =  (int) floor(oit->getColor() (GREEN) * std::abs(k) + AMB_LIGHT);
				if (g < COLOR_MAX)
					pg = g;
				else 
					pg = COLOR_MAX;
						
				int b = (int) floor(oit->getColor() (BLUE) * std::abs(k) + AMB_LIGHT);
				if (b < COLOR_MAX)
					pb = b;
				else 
					pb = COLOR_MAX;
						
				ublas::vector<float> fc = ray.getPixel();	
				fc (Z) = -c.getFarClip();
				float disN = norm_2(S - ray.getPixel());
				float disF = norm_2(fc - S);
	
				int pdepth = (int)  COLOR_MAX - std::min( float (COLOR_MAX), (COLOR_MAX)*
					(disN/(c.getFarClip() - c.getNearClip())));	
						
				img.setPixelRed(ray.getX(), ray.getY(), pr);
				img.setPixelGreen(ray.getX(), ray.getY(), pg);
				img.setPixelBlue(ray.getX(), ray.getY(), pb);
				img.setPixelDepth(ray.getX(), ray.getY(), pdepth);
				//std::cout<<ray.getX()<<","<<ray.getY()<<": "<<v<<" r^2:"<<oit->getRadiusSquared()
				//			<<" c^2:"<<oit->getDistanceToPRPSquared()<<" d:"<<d<<" "
				//			<<pr<<","<<pb<<","<<pg<<std::endl;
			}	//end if in view frustrum		
		}	//end if intersection possible	
	}	//end for each object/sphere
} 

