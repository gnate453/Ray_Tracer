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
			int i;
			while (parse>>i) {
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
			ublas::vector<float> vup (VECTOR_3DH);

			parse>>name;
			parse>>prp (X);
			parse>>prp (Y);
			parse>>prp (Z);
			prp (W) = 1.0;
			
			parse>>vpn (X);
			parse>>vpn (Y);
			parse>>vpn (Z);
			vpn (W) = 1.0;			

			parse>>vup (X);
			parse>>vup (Y);
			parse>>vup (Z);
			vup (W) = 1.0;
		
			float near, far;
			parse>>near;
			parse>>far;
			//TODO: parse new camera information.  Construct camera object appropriately.
			Camera camera(name, prp, vpn, vup, near, far);
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
			ublas::vector<float> n = (1/norm_2(cit->getVPN())) * cit->getVPN() ;
			ublas::vector<float> u = (1/norm_2(crossProductVectors(cit->getVUP(), n))) 
										* crossProductVectors(cit->getVUP(), n);
			ublas::vector<float> v = crossProductVectors(n, u);

			Image img(sit->getName()+"_"+cit->getName(), sit->getWidth(), sit->getHeight());
			for (int i = 0; i < sit->getHeight(); ++i)
			{
				
				for (int j = 0; j < sit->getWidth(); ++j)
				{
					float x = (2.0/(sit->getWidth()-1))*j - 1;
					float y = (2.0/(sit->getHeight()-1))*i - 1;
					ublas::vector<float> pixelWorldCoord (VECTOR_3DH);
					pixelWorldCoord (X) = cit->getPRP() (X) + (-cit->getNearClip() * cit->getVPN()(X)) + (x * u(X)) + (y * v(X));
					pixelWorldCoord (Y) = cit->getPRP() (Y) + (-cit->getNearClip() * cit->getVPN()(Y)) + (x * u(Y)) + (y * v(Y));
					pixelWorldCoord (Z) = cit->getPRP() (Z) + (-cit->getNearClip() * cit->getVPN()(Z)) + (x * u(Z)) + (y * v(Z));
					pixelWorldCoord (W) = 1.0;
					Ray tmpR(cit->getPRP(), pixelWorldCoord, j, i);
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
		//std::cout<<"PRP: "<<ray.getPRP()<<std::endl;
		//std::cout<<"VPN: "<<c.getVPN()<<std::endl;
		//std::cout<<"pixel coord: "<<ray.getPixelWorldCoord()<<std::endl;
		//std::cout<<"U: "<<ray.unitVector()<<std::endl;
		//c squared (sphere origin -camera prp)
		
		//std::cout<<" r: "<<s->getRadius()<<" c: "<<s->getDistanceToPixel(ray.getPixelWorldCoord())<<" v: "<<v<<std::endl;
		v = inner_prod(s->getOrigin() - ray.getPRP(), U);
		csqd = s->getDistanceToPixel(ray.getPRP()) *  s->getDistanceToPixel(ray.getPRP());
		dsqd = s->getRadiusSquared() - (csqd - (v * v));

		//std::cout<<" r^2: "<<s->getRadiusSquared()<<" c^2: "<<csqd<<" v^2: "<<v*v<<std::endl;
		if (dsqd >= 0) {
			float d = std::sqrt(dsqd);
			//std::cout<<"  d^2: "<<dsqd<<" d: "<<d<<"\n"<<std::endl;
			//normal to sphere
			ublas::vector<float> S = ray.paraPos(v-d);
			//if ( (norm_2(ray.paraPos(c.getNearClip())) < norm_2(S)) 
			//	&& (norm_2(S) < norm_2(ray.paraPos(c.getFarClip()))) ) { 
				
				ublas::vector<float> N = S - s->getOrigin();
				N = (1/norm_2(N)) * N; 
				
				float fr = 0;
				float fg = 0;
				float fb = 0;
				int sr = 0;
				int sg = 0;
				int sb = 0;
				int dr = 0;
				int dg = 0;
				int db = 0;
				
				for (std::list<Light>::iterator	l = lights.begin(); l != lights.end(); ++l) {
					
					//std::cout<<"Light: "<<l->getColor()<<std::endl;
					//std::cout<<"Object: "<<s->getColor().getSpecularProperties()<<std::endl;
					//std::cout<<"cos(theta): "<<dotProductVectors(N, l->getUnitVector())<<std::endl;
					
					//Specular lighting.
					double cosPhi = inner_prod(ray.unitVector(), ((2 * inner_prod(l->getUnitVector(), N) * N) - l->getUnitVector()));
					double phong = pow(cosPhi, s->getColor().getSpecularAlpha());
					
					//std::cout<<"cos(Phi): "<<cosPhi<<std::endl;
					//std::cout<<"Phong: "<<phong<<std::endl;

					fr = 0; 
					fg = 0;
					fb = 0;
					if (cosPhi <= 0) {
						fr = l->getRed() * s->getColor().getSpecularRed() * phong;
						fg = l->getGreen() * s->getColor().getSpecularGreen() * phong;
						fb = l->getBlue() * s->getColor().getSpecularBlue() * phong;
						
						if (fr > 0)	
							sr += (int) fr;
						if (fg > 0)
							sg += (int) fg;
						if (fb > 0)
							sb += (int) fb;
					}
					
					//Diffuse Lighting
					fr = 0; 
					fg = 0;
					fb = 0;
					float cosTheta = inner_prod(N, l->getUnitVector()); 
					if (cosTheta >= 0) {
						fr =  DIFFUSE_FACT * s->getColor().getDiffuseRed() 
										* std::abs(cosTheta);
						fg =  DIFFUSE_FACT * s->getColor().getDiffuseGreen() 
										* std::abs(cosTheta);
						fb =  DIFFUSE_FACT * s->getColor().getDiffuseBlue() 
										* std::abs(cosTheta);
						dr += (int) floor(((float) l->getRed()) * fr);
						dg += (int) floor(((float) l->getGreen()) * fg);
						db += (int) floor(((float) l->getBlue()) * fb);
					}
				}			
				
				//std::cout<<"Object: "<<s->getName()<<" "<<ray.getScreenX()<<" "<<ray.getScreenY()<<std::endl;
				//std::cout<<"red: "<<sr<<" "<<dr<<std::endl;
				//std::cout<<"green: "<<sg<<" "<<dg<<std::endl;
				//std::cout<<"blue: "<<sb<<" "<<db<<std::endl;


				//float k = inner_prod(ray.unitVector() , N);
				//std::cout<<"  Q: "<<Q<<" k: "<<k<<std::endl;


				//Ambient Lighting.
				int ar = s->getColor().getAmbientRed() * AMB_LIGHT;
				int ag = s->getColor().getAmbientGreen() * AMB_LIGHT;
				int ab = s->getColor().getAmbientBlue() * AMB_LIGHT;
					
				//add lighting.				
				pr = sr + dr + ar;
				if ( std::abs(pr) > COLOR_MAX)
					pr = COLOR_MAX;

				pg = sg + dg + ag;
				if ( std::abs(pg) > COLOR_MAX)
					pg = COLOR_MAX;

				pb = sb + db + ab;
				if ( std::abs(pb) > COLOR_MAX)
					pb = COLOR_MAX;

				//depth TODO: needs fix?	
				ublas::vector<float> fc = ray.getPixelWorldCoord();	
				fc (Z) = -c.getFarClip();
				float disN = norm_2(S - ray.getPixelWorldCoord());
				float disF = norm_2(fc - S);

				int pdepth = (int)  COLOR_MAX - std::min( float (COLOR_MAX), (COLOR_MAX)*
					(disN/(c.getFarClip() - c.getNearClip())));

				img.setPixelRed(ray.getScreenX(), ray.getScreenY(), pr);
				img.setPixelGreen(ray.getScreenX(), ray.getScreenY(), pg);
				img.setPixelBlue(ray.getScreenX(), ray.getScreenY(),  pb);
				img.setPixelDepth(ray.getScreenX(), ray.getScreenY(), pdepth);

			//}	//end if in view frustrum		
		}	//end if intersection possible
	}	//end for each sphere
} 


void intersectRayWithPolygons(Ray ray, std::list<Polygon> polygons, std::list<Light> lights, Camera c, Image &img) {
	
}
