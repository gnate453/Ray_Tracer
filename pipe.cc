#include "pipe.h"

World worldFromString(std::string d) {
	
	World newWorld;
	std::list<Face> facesInGroup;
	std::string currentGroupName = "Default";
	std::string curMaterial;
	std::stringstream data(d);
	std::string tmp;
	bool facesAdded = false;
	bool groupNamed = false;

	while (std::getline(data, tmp))
	{
		//vertex
		if (*(tmp.begin()) == 'v') {
			std::stringstream parse(tmp.substr(LENGTH_1, tmp.length()));
			ublas::vector<float> v (VECTOR_3DH);
			parse>>v(X);
			parse>>v(Y);
			parse>>v(Z);

			float n;
			if (parse>>n)
				v(W) = n;
			else
				v(W) = 1.0;			

			//std::cout<<"v: "<<v<<std::endl;			

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
				a = points.front() - 1;
				points.pop();
				b = points.front() - 1;
				points.pop();
				c = points.front() - 1;
				points.pop();
				Face cFace(newWorld.getVertices()[a], 
							newWorld.getVertices()[b], 
							newWorld.getVertices()[c]);
				facesInGroup.push_back(cFace);
				//std::cout<<"f: "<<newWorld.getVertices()[a]<< 
				//			newWorld.getVertices()[b]<< 
				//			newWorld.getVertices()[c]<<std::endl;			
			}
			else {
				int a = points.front() - 1;
				points.pop();
				while (points.size() >= (VERTEX_PER_FACE - 1)) {
					int b, c;
					b = points.front() - 1;
					points.pop();
					c = points.front() - 1;
					Face cFace(newWorld.getVertices()[a], 
								newWorld.getVertices()[b], 
								newWorld.getVertices()[c]);
					facesInGroup.push_back(cFace);
					//std::cout<<"f: "<<newWorld.getVertices()[a]<< 
					//			newWorld.getVertices()[b]<< 
					//			newWorld.getVertices()[c]<<std::endl;			
				}
			}
			
			facesAdded = true;
		}
		//group
		else if (*(tmp.begin()) == 'g') {
			std::stringstream parse(tmp.substr(LENGTH_1, tmp.length()));

			if ((!facesAdded && groupNamed) || (!facesAdded && !groupNamed)){
				//second group, default or no faces in last group.
				//keep list of faces
				//get new group name
				parse>>currentGroupName;
			}
			else if ((facesAdded && !groupNamed) || (facesAdded && groupNamed)) {		
				Polygon p(currentGroupName, (newWorld.getMaterials())[curMaterial], facesInGroup);
				newWorld.addPolygon(p);
				parse>>currentGroupName;
				facesInGroup.clear();
			}

			groupNamed = true;			
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

	if ( (facesAdded && !groupNamed) || (facesAdded && groupNamed) ) {
		if (!facesInGroup.empty()) {
			Polygon p(currentGroupName, (newWorld.getMaterials())[curMaterial], facesInGroup);
			newWorld.addPolygon(p);
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
		int imgCount = 0;
		for (std::list<Camera>::iterator cit = cameras.begin();
				cit != cameras.end(); ++cit)
		{	
			ublas::vector<float> camU = cit->getHorizontalVector();
			ublas::vector<float> camV = cit->getVerticalVector();

			int count = 0;
			int percent = 0;
			Image* img = new Image(sit->getName()+"_"+cit->getName(), sit->getWidth(), sit->getHeight());
			for (int i = 0; i < sit->getHeight(); ++i)
			{
				
				for (int j = 0; j < sit->getWidth(); ++j)
				{
					float x = (2.0/(sit->getWidth()-1))*j - 1;
					float y = (2.0/(sit->getHeight()-1))*i - 1;
					
					ublas::vector<float> r(VECTOR_3DH);
					ublas::vector<float> pixelWorldCoord (VECTOR_3DH);
					r (X) = (-cit->getNearClip() * cit->getVPN()(X)) + (x * camU(X)) + (y * camV(X));
					r (Y) = (-cit->getNearClip() * cit->getVPN()(Y)) + (x * camU(Y)) + (y * camV(Y));
					r (Z) = (-cit->getNearClip() * cit->getVPN()(Z)) + (x * camU(Z)) + (y * camV(Z));
					r (W) = 1.0;
					pixelWorldCoord (X) = cit->getPRP() (X) + r (X);
					pixelWorldCoord (Y) = cit->getPRP() (Y) + r (Y);
					pixelWorldCoord (Z) = cit->getPRP() (Z) + r (Z);
					pixelWorldCoord (W) = 1.0;
					Ray tmpR(cit->getPRP(), r, pixelWorldCoord, j, i, x, y);

					intersectRayWithSpheres(tmpR, spheres,	w.getLights(), *cit, *img);
					intersectRayWithPolygons(tmpR, polygons, w.getLights(), *cit, *img);
					count += 1;
					if (count % (sit->getWidth()*sit->getHeight() / 10) == 0){
						percent += 10;
						std::cout<<"Image "<<imgCount<<" progress: "<<percent<<"%"<<std::endl;
					}
				}
				//std::cout<<std::endl;
			}	
					
			imgs.push_back(*img);
			imgCount += 1;
		}//end for each camera
	}//end for each scene


	outputImages(imgs);
}

void intersectRayWithSpheres(Ray ray, std::list<Sphere> spheres, std::list<Light> lights, Camera c, Image &img) {

	ublas::vector<float> U = ray.unitVector();
	//std::cout<<ray.getX()<<","<<ray.getY()<<" U: "<<U<<std::endl;
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
			if ( ray.paraPos(c.getNearClip())(Z) > S (Z) 
				&& S (Z) > ray.paraPos(c.getFarClip()) (Z) ) { 
				
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
				int pr = 0;
				int pg = 0;
				int pb = 0;
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
						fr = SPECULAR_REDUCT * l->getRed() * s->getColor().getSpecularRed() * phong;
						fg = SPECULAR_REDUCT * l->getGreen() * s->getColor().getSpecularGreen() * phong;
						fb = SPECULAR_REDUCT * l->getBlue() * s->getColor().getSpecularBlue() * phong;
						
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

			}	//end if in view frustrum		
		}	//end if intersection possible
	}	//end for each sphere
} 


void intersectRayWithPolygons(Ray ray, std::list<Polygon> polygons, std::list<Light> lights, Camera c, Image &img) {
	ublas::vector<float> U = ray.unitVector();
	ublas::vector<float> d = (ray.getAlpha() * c.getHorizontalVector()) 
								+ (ray.getBeta() * c.getVerticalVector())
								- (c.getNearClip() * c.getVPN());
	//std::cout<<"U: "<<U<<" d:"<<d<<std::endl;
	for	(std::list<Polygon>::iterator p = polygons.begin(); p != polygons.end(); ++p) {
		int hitCount = 0;
		//Face *closestFace;
		float disToClosestFace = c.getFarClip();
		std::list<Face> facesInP = p->getFaces();
		for (std::list<Face>::iterator f = facesInP.begin(); f != facesInP.end(); ++f) {
			if (inner_prod(ray.unitVector(), f->getNormal()) != 0) {
				float tstar = (1/inner_prod(ray.unitVector(), f->getNormal())) 
								*  (-inner_prod(ray.getPRP(), f->getNormal() + d));
				
				//princeton.edu cs426 fall 00 lectures on raycasting
				ublas::vector<float> S = ray.paraPos(tstar);
				ublas::vector<float> V1 = f->getVertex(P_ONE) - S;
				ublas::vector<float> V2 = f->getVertex(P_TWO) - S;
				ublas::vector<float> N1 = crossProductVectors(V2, V1);
				N1 = (1/norm_2(N1)) * N1;
				float d1 = inner_prod(-ray.getPRP(), N1);
				bool on1 = ( (inner_prod(S, N1) + d1) < 0);				

				ublas::vector<float> V3 = f->getVertex(P_TWO) - S;
				ublas::vector<float> V4 = f->getVertex(P_THREE) - S;
				ublas::vector<float> N2 = crossProductVectors(V4, V3);
				N2 = (1/norm_2(N2)) * N2;
				float d2 = inner_prod(-ray.getPRP(), N2);				
				bool on2 = ( (inner_prod(S, N2) + d2) < 0 );				

				ublas::vector<float> V5 = f->getVertex(P_THREE) - S;
				ublas::vector<float> V6 = f->getVertex(P_ONE) - S;
				ublas::vector<float> N3 = crossProductVectors(V6, V5);
				N3 = (1/norm_2(N3)) * N3;
				float d3 = inner_prod(-ray.getPRP(), N3);
				bool on3 = ( (inner_prod(S, N3) + d3) < 0 );

				if (on1 && on2 && on3) {
					float fr = 0;
					float fg = 0;
					float fb = 0;
					int sr = 0;
					int sg = 0;
					int sb = 0;
					int dr = 0;
					int dg = 0;
					int db = 0;
					int pr = 0;
					int pg = 0;
					int pb = 0;
					for (std::list<Light>::iterator	l = lights.begin(); l != lights.end(); ++l) {
					
						//Specular lighting.
						double cosPhi = inner_prod(ray.unitVector(), ((2 * inner_prod(l->getUnitVector(), 
							f->getNormal()) * f->getNormal()) - l->getUnitVector()));
						double phong = pow(cosPhi, p->getColor().getSpecularAlpha());
				
						fr = 0; 
						fg = 0;
						fb = 0;
						if (cosPhi <= 0) {
							fr = l->getRed() * p->getColor().getSpecularRed() * phong;
							fg = l->getGreen() * p->getColor().getSpecularGreen() * phong;
							fb = l->getBlue() * p->getColor().getSpecularBlue() * phong;
							
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
						float cosTheta = inner_prod(f->getNormal(), l->getUnitVector()); 
						if (cosTheta >= 0) {
							fr =  DIFFUSE_FACT * p->getColor().getDiffuseRed() 
											* std::abs(cosTheta);
							fg =  DIFFUSE_FACT * p->getColor().getDiffuseGreen() 
											* std::abs(cosTheta);
							fb =  DIFFUSE_FACT * p->getColor().getDiffuseBlue() 
											* std::abs(cosTheta);
							dr += (int) floor(((float) l->getRed()) * fr);
							dg += (int) floor(((float) l->getGreen()) * fg);
							db += (int) floor(((float) l->getBlue()) * fb);
						}
					}			
						
					//Ambient Lighting.
					int ar = p->getColor().getAmbientRed() * AMB_LIGHT;
					int ag = p->getColor().getAmbientGreen() * AMB_LIGHT;
					int ab = p->getColor().getAmbientBlue() * AMB_LIGHT;
							
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
				}			
			}
	
		} //end for each face
	} //end for each polygon
}



						


