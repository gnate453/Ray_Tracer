#include "pipe.h"

World worldFromString(std::string d) {
	
	World newWorld;
	std::list<Face> facesInGroup;
	std::string currentGroupName = "Default";
	std::string curMaterial;
	std::stringstream data(d);
	std::string tmp;

	while (std::getline(data, tmp))
	{
		//vertex
		if (*(tmp.begin()) == 'v') {
			std::stringstream parse(tmp.substr(LENGTH_1, tmp.length()));
			ublas::vector<float> v (VECTOR_3DH);
			parse>>v(X);
			parse>>v(Y);
			parse>>v(Z);

			if (!parse.str().empty())
				parse>>v(W);
			else
				v(W) = 1.0;			

			newWorld.addVertex(v);
		}
		//face
		else if (*(tmp.begin()) == 'f') {
			std::stringstream parse(tmp.substr(LENGTH_1, tmp.length()));
			
			std::queue<int> points;
			while (!parse.str().empty()) {
				char v;
				parse>>v;
				int i;
				parse>>i;
				points.push(i);
			}
			
			if (points.size() == VERTEX_PER_FACE) {
				int a, b, c;
				a = points.front();
				points.pop();
				b = points.front();
				points.pop();
				c = points.front();
				points.pop();
				Face cFace(newWorld.getVertices()[a], 
							newWorld.getVertices()[b], 
							newWorld.getVertices()[c]);
				facesInGroup.push_back(cFace);
			}
			else {
				int a = points.front();
				points.pop();
				while (points.size() >= (VERTEX_PER_FACE - 1)) {
					int b, c;
					b = points.front();
					points.pop();
					c = points.front();
					Face cFace(newWorld.getVertices()[a], 
								newWorld.getVertices()[b], 
								newWorld.getVertices()[c]);
					facesInGroup.push_back(cFace);
				}
			}
			
		}
		//group
		else if (*(tmp.begin()) == 'g') {
			std::stringstream parse(tmp.substr(LENGTH_1, tmp.length()));

			if (!facesInGroup.empty())
			{
				Polygon p(currentGroupName, (newWorld.getMaterials())[curMaterial], facesInGroup);
				newWorld.addPolygon(p);
			}
			
			parse>>currentGroupName;
			facesInGroup = std::list<Face>();
		}
		//sphere
		else if (*(tmp.begin()) == 's') {
			std::stringstream parse(tmp.substr(LENGTH_1, tmp.length()));
			
			std::string name;
			ublas::vector<float> w (VECTOR_3DH);
			
			parse>>name;
			parse>>w (X);
			parse>>w (Y);
			parse>>w (Z);
			w (W) = 1.0;
			float r;
			parse>>r;
			//std::cout<<name<<" "<<w<<" "<<r<<std::endl;
			//TODO: change sphere to have material instead of color.
			Sphere sphere(name, newWorld.getMaterials()[curMaterial], w, r);
			newWorld.addSphere(sphere);
		}
		//new material
		else if (*(tmp.begin()) == 'n') { 
			std::stringstream parse(tmp.substr(LENGTH_6, tmp.length()));
		
			//std::cout<<parse.str()<<std::endl;	
			std::string name, prefix;
			ublas::vector<float> ka (VECTOR_C);
			ublas::vector<float> kd (VECTOR_C);
			ublas::vector<float> ks (VECTOR_C);

			parse>>name;

			std::getline(data, tmp);
			std::stringstream a(tmp.substr(LENGTH_2, tmp.length()));
			a>>ka(RED);
			a>>ka(GREEN);
			a>>ka(BLUE);
			ka(ALPHA) = 0.0;
			std::getline(data, tmp);
			std::stringstream d(tmp.substr(LENGTH_2, tmp.length()));
			d>>kd(RED);
			d>>kd(GREEN);
			d>>kd(BLUE);
			kd(ALPHA) = 0.0;
			std::getline(data, tmp);
			std::stringstream s(tmp.substr(LENGTH_2, tmp.length()));
			s>>ks(RED);
			s>>ks(GREEN);
			s>>ks(BLUE);
			std::getline(data, tmp);
			std::stringstream ph(tmp.substr(LENGTH_2, tmp.length()));
			ph>>ks(ALPHA);
			
			//std::cout<<name<<" "<<ka<<" "<<kd<<" "<<ks<<std::endl; 			

			Material mtl(name, ka, kd, ks);
			newWorld.addMaterial(name, mtl);
		}
		//use material
		else if (*(tmp.begin()) == 'u') {	
			std::stringstream parse(tmp.substr(LENGTH_6, tmp.length()));
			//std::cout<<"Before use: "<<curMaterial<<std::endl;	
			parse>>curMaterial;
			//std::cout<<"After use: "<<curMaterial<<std::endl;	
		}
		//camera
		else if (*(tmp.begin()) == 'c')
		{	
			std::stringstream parse(tmp.substr(LENGTH_1, tmp.length()));
			
			std::string name;
			ublas::vector<float> prp (VECTOR_3DH);
			ublas::vector<float> vpn (VECTOR_3DH);

			parse>>name;
			parse>>prp (X);
			parse>>prp (Y);
			parse>>prp (Z);
			prp (W) = 1.0;
			
			parse>>vpn (X);
			parse>>vpn (Y);
			parse>>vpn (Z);
			vpn (W) = 1.0;			

			float near, far;
			parse>>near;
			parse>>far;
			//put the vpn at z= -near;
			vpn (Z) = -near; 
			//TODO: add new camera information to Camera Class.
			//TODO: parse new camera information.  Construct camera object appropriately.
			Camera camera(name, prp, vpn, near, far);
			newWorld.addCamera(camera);
		}
		//scene
		else if (*(tmp.begin()) == 'r')
		{	
			//TODO: DOUBLE CHECK this is still all that is needed for a scene.
			std::stringstream parse(tmp.substr(LENGTH_1, tmp.length()));

			std::string name;
			int width, height, depth;
			
			parse>>name;
			parse>>width;
			parse>>height;
			parse>>depth;

			Scene scene(name, width, height, depth);
			newWorld.addScene(scene);
		}
		else if (*(tmp.begin()) == 'l')
		{	
			std::stringstream parse(tmp.substr(LENGTH_1, tmp.length()));

			ublas::vector<float> dir (VECTOR_3DH);
			ublas::vector<int> c (VECTOR_C);

			parse>>dir(X);
			parse>>dir(Y);
			parse>>dir(Z);
			parse>>dir(W);
			parse>>c(RED);
			parse>>c(GREEN);
			parse>>c(BLUE);
			c(ALPHA) = 1;		
		
			Light l(dir, c);
			
			newWorld.addLight(l);
		}
	}

	//std::cout<<"parse 2 done"<<std::endl;

	return newWorld;
}

void castRays(World w) {

	std::list<Sphere> spheres = w.getSpheres();
	std::list<Polygon> polygons = w.getPolygons();
	std::list<Camera> cameras = w.getCameras();
	std::list<Scene> scenes = w.getScenes();
	std::list<Image> imgs;
	//float rsqd, csqd;
	
	for (std::list<Scene>::iterator sit = scenes.begin();
			sit != scenes.end(); ++sit)
	{
		for (std::list<Camera>::iterator cit = cameras.begin();
				cit != cameras.end(); ++cit)
		{
			//get r^2 and c^2 for each sphere.
			//for (std::list<Sphere>::iterator oit = spheres.begin();
			//		oit != spheres.end(); ++oit)
			//{
			//	std::cout<<oit->getName()<<": "<<oit->getOrigin()<<" "<<oit->getRadius()<<std::endl;	
			//}//end for each sphere

			Image img(sit->getName()+"_"+cit->getName(), sit->getWidth(), sit->getHeight());
			for (int i = 0; i < sit->getHeight(); ++i)
			{
				float y = (2.0/(sit->getHeight()-1))*i - 1 ;
				for (int j = 0; j < sit->getWidth(); ++j)
				{
					Ray tmpR(cit->getPRP());
					float x = (2.0/(sit->getWidth()-1))*j - 1;
					tmpR.setScreenX(j);
					tmpR.setScreenY(i);
					//std::cout<<x<<","<<y<<" ";
					//std::cout<<j<<","<<i<<" ";
					ublas::vector<float> rv (VECTOR_3DH);
					rv (X) = x;
					rv (Y) = y;
					rv (Z) = -cit->getNearClip();
					rv (W) = 1.0;
					tmpR.setPixelWorldCoord(rv);
					intersectRayWithSpheres(tmpR, spheres,	w.getLights(), *cit, img);
					intersectRayWithPolygons(tmpR, polygons, w.getLights(), *cit, img);
				}
				//std::cout<<std::endl;
			}	
					
			imgs.push_back(img);

		}//end for each camera
	}//end for each scene


	outputImages(imgs);
}

void intersectRayWithSpheres(Ray ray, std::list<Sphere> spheres, std::list<Light> lights, Camera c, Image &img) {

	ublas::vector<float> U = ray.unitVector();
	//std::cout<<ray.getX()<<","<<ray.getY()<<" U: "<<U<<std::endl;
	int pr = 0;
	int pg = 0;
	int pb = 0;
	for (std::list<Sphere>::iterator s = spheres.begin(); s != spheres.end(); ++s)
	{
		float v, csqd, dsqd;
		//std::cout<<"Object: "<<s->getName()<<" "<<ray.getScreenX()<<" "<<ray.getScreenY()<<std::endl;
		//c squared (sphere origin -camera prp)
		if (ray.getPixelWorldCoord() (Z) > s->getOrigin() (Z)) {
		//	std::cout<<"sphere past near clip "<<std::endl;
			csqd = s->getDistanceToVPN(ray.getPixelWorldCoord()) *  s->getDistanceToVPN(ray.getPixelWorldCoord());
			v = inner_prod(subtractVectors(s->getOrigin(), ray.getPixelWorldCoord()), U);
			dsqd = s->getRadiusSquared() - ( csqd - (v * v) );
			//std::cout<<dsqd<<std::endl;
		}
		else {
			//	std::cout<<"sphere before near clip "<<std::endl;
			csqd = s->getDistanceToVPN(ray.getPixelWorldCoord()) *  s->getDistanceToVPN(ray.getPixelWorldCoord());
			v = inner_prod(subtractVectors(ray.getPixelWorldCoord(),  s->getOrigin()), U);
			dsqd = ( csqd + (v * v) ) - s->getRadiusSquared();
			//std::cout<<dsqd<<std::endl;
		}
		//std::cout<<" r^2: "<<s->getRadiusSquared()<<" c^2: "<<csqd<<" v^2: "<<v*v<<std::endl;
		if (dsqd >= 0)
		{
			float d = std::sqrt(dsqd);
			//std::cout<<"  d^2: "<<dsqd<<" d: "<<d<<"\n"<<std::endl;
			//normal to sphere
			ublas::vector<float> S = ray.paraPos(v-d);
			if ( S(Z) < -c.getNearClip()  && S(Z) > -c.getFarClip()) { 
				//std::cout<<"Object: "<<oit->getName()<<" S(Z): "<<S(Z)<<std::endl;
				//std::cout<<"near: "<<-cit->getNearClip()<<" far: "<<-cit->getFarClip()<<std::endl;
				ublas::vector<float> N = subtractVectors(s->getOrigin(), S);
				N = (1/norm_2(N)) * N; 
				

				float fr = 0;
				float fg = 0;
				float fb = 0;
				int dr = 0;
				int dg = 0;
				int db = 0;
				for (std::list<Light>::iterator	l = lights.begin(); l != lights.end(); ++l) {
					//std::cout<<"Light: "<<l->getColor()<<std::endl;
					//std::cout<<"Object: "<<s->getColor().getDiffuseProperties()<<std::endl;
					//std::cout<<"cos(theta): "<<inner_prod(N, l->getUnitVector())<<std::endl;
					if (inner_prod(N, l->getUnitVector()) < 0) {
						fr =  DIFFUSE_FACT * s->getColor().getDiffuseRed() 
										* std::abs(inner_prod(N, l->getUnitVector()));
						fg =  DIFFUSE_FACT * s->getColor().getDiffuseGreen() 
										* std::abs(inner_prod(N, l->getUnitVector()));
						fb =  DIFFUSE_FACT * s->getColor().getDiffuseBlue() 
										* std::abs(inner_prod(N, l->getUnitVector()));
						dr += (int) floor(((float) l->getRed()) * fr);
						dg += (int) floor(((float) l->getGreen()) * fg);
						db += (int) floor(((float) l->getBlue()) * fb);
						}
				}			

				
				//std::cout<<"Object: "<<s->getName()<<" "<<ray.getScreenX()<<" "<<ray.getScreenY()<<std::endl;
				//std::cout<<"red: "<<dr<<std::endl;
				//std::cout<<"green: "<<dg<<std::endl;
				//std::cout<<"blue: "<<db<<std::endl;

				//float k = inner_prod(ray.unitVector() , N);
				//std::cout<<"  Q: "<<Q<<" k: "<<k<<std::endl;

				int ar = s->getColor().getAmbientRed() * AMB_LIGHT;
				int ag = s->getColor().getAmbientGreen() * AMB_LIGHT;
				int ab = s->getColor().getAmbientBlue() * AMB_LIGHT;
					
				ublas::vector<float> fc = ray.getPixelWorldCoord();	
				fc (Z) = -c.getFarClip();
				float disN = norm_2(subtractVectors(S, ray.getPixelWorldCoord()));
				float disF = norm_2(subtractVectors(fc, S));

				int pdepth = (int)  COLOR_MAX - std::min( float (COLOR_MAX), (COLOR_MAX)*
					(disN/(c.getFarClip() - c.getNearClip())));	
				
				
				pr = dr + ar;
				if (pr > COLOR_MAX)
					pr = COLOR_MAX;

				pg = dg + ag;
				if (pg >COLOR_MAX)
					pg = COLOR_MAX;

				pb = db + ab;
				if (pb > COLOR_MAX)
					pb = COLOR_MAX;
		
				img.setPixelRed(ray.getScreenX(), ray.getScreenY(), pr);
				img.setPixelGreen(ray.getScreenX(), ray.getScreenY(), pg);
				img.setPixelBlue(ray.getScreenX(), ray.getScreenY(),  pb);
				img.setPixelDepth(ray.getScreenX(), ray.getScreenY(), pdepth);
				//std::cout<<ray.getX()<<","<<ray.getY()<<": "<<v<<" r^2:"<<oit->getRadiusSquared()
				//			<<" c^2:"<<oit->getDistanceToPRPSquared()<<" d:"<<d<<" "
				//			<<pr<<","<<pb<<","<<pg<<std::endl;
			}	//end if in view frustrum		
		}	//end if intersection possible
	}	//end for each sphere
} 


void intersectRayWithPolygons(Ray ray, std::list<Polygon> polygons, std::list<Light> lights, Camera c, Image &img) {
	
}
