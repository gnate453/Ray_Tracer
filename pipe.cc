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
			ublas::vector<float> v (VECTOR_3D);
			parse>>v(X);
			parse>>v(Y);
			parse>>v(Z);

			float n;
			if (parse>>n) {
				v(X) = v(X)/n;
				v(Y) = v(Y)/n;
				v(Z) = v(Z)/n;
			}

			std::cout<<"v: "<<v<<std::endl;			

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
				std::cout<<"f: "<<newWorld.getVertices()[a]<< 
							newWorld.getVertices()[b]<< 
							newWorld.getVertices()[c]<<std::endl;			
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
					std::cout<<"f: "<<newWorld.getVertices()[a]<< 
								newWorld.getVertices()[b]<< 
								newWorld.getVertices()[c]<<std::endl;			
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
			ublas::vector<float> w (VECTOR_3D);
			
			parse>>name;
			parse>>w (X);
			parse>>w (Y);
			parse>>w (Z);
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
			ublas::vector<float> prp (VECTOR_3D);
			ublas::vector<float> vpn (VECTOR_3D);
			ublas::vector<float> vup (VECTOR_3D);

			parse>>name;
			parse>>prp (X);
			parse>>prp (Y);
			parse>>prp (Z);
			
			parse>>vpn (X);
			parse>>vpn (Y);
			parse>>vpn (Z);		

			parse>>vup (X);
			parse>>vup (Y);
			parse>>vup (Z);

		
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

			ublas::vector<float> dir (VECTOR_3D);
			ublas::vector<int> c (VECTOR_C);
			float w;

			parse>>dir(X);
			parse>>dir(Y);
			parse>>dir(Z);
			parse>>w;

			dir(X) = dir(X)/w;
			dir(Y) = dir(Y)/w;
			dir(Z) = dir(Z)/w; 
			
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

			//std::cout<<"U: "<<camU<<" V: "<<camV<<" PRP: "<<cit->getPRP()<<std::endl;

			int count = 0;
			int percent = 0;
			Image* img = new Image(sit->getName()+"_"+cit->getName(), sit->getWidth(), sit->getHeight());
			for (int i = 0; i < sit->getHeight(); ++i)
			{
				
				for (int j = 0; j < sit->getWidth(); ++j)
				{
					float x = (2.0/(sit->getWidth()-1))*j - 1;
					float y = (2.0/(sit->getHeight()-1))*i - 1;
					
					ublas::vector<float> r(VECTOR_3D);
					ublas::vector<float> pixelWorldCoord (VECTOR_3D);

					r (X) = (-cit->getNearClip() * cit->getVPN()(X)) + (x * camU(X)) + (y * camV(X));
					r (Y) = (-cit->getNearClip() * cit->getVPN()(Y)) + (x * camU(Y)) + (y * camV(Y));
					r (Z) = (-cit->getNearClip() * cit->getVPN()(Z)) + (x * camU(Z)) + (y * camV(Z));

					pixelWorldCoord (X) = cit->getPRP() (X) + r (X);
					pixelWorldCoord (Y) = cit->getPRP() (Y) + r (Y);
					pixelWorldCoord (Z) = cit->getPRP() (Z) + r (Z);
					
					//std::cout<<"Ray: "<<r<<" Pixel: "<<pixelWorldCoord<<std::endl;
		
					Ray tmpR(cit->getPRP(), r, pixelWorldCoord, j, i, x, y);
		
					if (!spheres.empty() && !polygons.empty()) {
						Intersection closestSphere = intersectRayWithSpheres(tmpR, spheres, *cit);
						Intersection closestPolygon = intersectRayWithPolygons(tmpR, polygons, *cit);
						Intersection* closest;
						if (closestSphere.getDepth() < closestPolygon.getDepth())
							closest = &closestSphere;
						else
							closest = &closestPolygon;

						if ( cit->getNearClip() < std::abs(closest->getDepth()) && std::abs(closest->getDepth()) < cit->getFarClip() ) {
							ublas::vector<int> color;
							color = calcPixelColor(tmpR, closest->getPoint(), closest->getSurfaceNormal(),
									 	closest->getSurfaceMaterial(), w.getLights(), cit->getNearClip(), cit->getFarClip());

							img->setPixelRed(j, i, color(RED));
							img->setPixelGreen(j, i, color(GREEN));
							img->setPixelBlue(j, i,  color(BLUE));
							img->setPixelDepth(j, i, color(ALPHA));
						}
					}
					else if (!spheres.empty() &&  polygons.empty()) {
						Intersection closest = intersectRayWithSpheres(tmpR, spheres, *cit);
						if ( cit->getNearClip() < std::abs(closest.getDepth()) && std::abs(closest.getDepth()) < cit->getFarClip() ) {
							ublas::vector<int> color;
							color = calcPixelColor(tmpR, closest.getPoint(), closest.getSurfaceNormal(),
									 	closest.getSurfaceMaterial(), w.getLights(), cit->getNearClip(), cit->getFarClip());

							img->setPixelRed(j, i, color(RED));
							img->setPixelGreen(j, i, color(GREEN));
							img->setPixelBlue(j, i,  color(BLUE));
							img->setPixelDepth(j, i, color(ALPHA));
						}
					}
					else if (spheres.empty() &&  !polygons.empty()) {
						Intersection closest = intersectRayWithPolygons(tmpR, polygons, *cit);
						if ( cit->getNearClip() < std::abs(closest.getDepth()) && std::abs(closest.getDepth()) < cit->getFarClip() ) {
							ublas::vector<int> color;
							color = calcPixelColor(tmpR, closest.getPoint(), closest.getSurfaceNormal(),
									 	closest.getSurfaceMaterial(), w.getLights(), cit->getNearClip(), cit->getFarClip());

							img->setPixelRed(j, i, color(RED));
							img->setPixelGreen(j, i, color(GREEN));
							img->setPixelBlue(j, i,  color(BLUE));
							img->setPixelDepth(j, i, color(ALPHA));
						}
					}


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

	std::cout<<"Writing image files..."<<std::endl;
	outputImages(imgs);
	std::cout<<"Done!"<<std::endl;
}

Intersection intersectRayWithSpheres(Ray ray, std::list<Sphere> spheres, Camera c) {

	ublas::vector<float> U = ray.unitVector();
	float closest = c.getFarClip();
	Intersection closestIntersection = Intersection(closest, ray.paraPos(closest), U, spheres.front().getColor());
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
		csqd = s->getDistanceToPixel(ray.getPRP()) * s->getDistanceToPixel(ray.getPRP());
		dsqd = s->getRadiusSquared() - (csqd - (v * v));

		//std::cout<<" r^2: "<<s->getRadiusSquared()<<" c^2: "<<csqd<<" v^2: "<<v*v<<std::endl;
		if (dsqd >= 0) {
			float d = std::sqrt(dsqd);
			//std::cout<<"  d^2: "<<dsqd<<" d: "<<d<<"\n"<<std::endl;
			//normal to sphere
			ublas::vector<float> S = ray.paraPos(v-d);
			ublas::vector<float> toNear = S - ray.paraPos(c.getNearClip());			
			ublas::vector<float> toFar = ray.paraPos(c.getFarClip()) - S;

			if ( c.getNearClip() < std::abs(v-d) && std::abs(v-d) < c.getFarClip() ) { 
				
				ublas::vector<float> N = S - s->getOrigin();
				N = (1/norm_2(N)) * N; 
				
				//record closest distance interestion with infromation to color pixel.
				if ((v-d) < closest) {
					closest = v-d;
					closestIntersection = Intersection( (v-d), S, N, s->getColor() );
				}	

			}	//end if in view frustrum		
		}	//end if intersection possible
	}	//end for each sphere

	return closestIntersection;
} 


Intersection intersectRayWithPolygons(Ray ray, std::list<Polygon> polygons, Camera c) {
	ublas::vector<float> U = ray.unitVector();
	
	float closest = c.getFarClip();
	Intersection closestIntersection = Intersection(closest, ray.paraPos(closest), U, polygons.front().getColor());

	for	(std::list<Polygon>::iterator p = polygons.begin(); p != polygons.end(); ++p) {
		int hitCount = 0;
		//Face *closestFace;
		float disToClosestFace = c.getFarClip();
		std::list<Face> facesInP = p->getFaces();
		for (std::list<Face>::iterator f = facesInP.begin(); f != facesInP.end(); ++f) {
			ublas::vector<float> e1 = f->getVertex(P_TWO) - f->getVertex(P_ONE);
			ublas::vector<float> e2 = f->getVertex(P_THREE) - f->getVertex(P_ONE);
			
			ublas::vector<float> prod = crossProductVectors(U, e2);

			float det = inner_prod(e1, prod);
			
			if (det < -EPSILON || EPSILON < det) {
				det = 1.0/det;
				
				ublas::vector<float> T = f->getVertex(P_ONE) - ray.getPRP();
				float u = inner_prod(T, prod) * det;

				if ( 0.0 <= u && u <= 1.0) {
					ublas::vector<float> Q = crossProductVectors(T, e1);

					float v = inner_prod(U, Q) * det;
					if ( 0.0 <= v && u+v <= 1.0) {
						float tstar = inner_prod(e2, Q) * det;
						if (tstar > EPSILON) {
							if (tstar < closest){
								closest = tstar;
								closestIntersection = Intersection(closest, ray.paraPos(tstar), f->getNormal(), p->getColor());
							}
						} //end if intersection 					
					} //end if v in range
				} //end if u in range 
			} //end if det not close to zero 			
		} //end for each face
	} //end for each polygon

	return closestIntersection;
}

ublas::vector<float> calcPixelColor(Ray ray, ublas::vector<float> point, ublas::vector<float> surfaceNormal, Material surfaceMaterial, std::list<Light> lights, float near, float far) {

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
		double cosPhi = inner_prod(ray.unitVector(), ((2 * inner_prod(l->getUnitVector(),
									 surfaceNormal) * surfaceNormal) - l->getUnitVector()));
		double phong = pow(cosPhi, surfaceMaterial.getSpecularAlpha());

		//std::cout<<"cos(Phi): "<<cosPhi<<std::endl;
		//std::cout<<"Phong: "<<phong<<std::endl;

		fr = 0; 
		fg = 0;
		fb = 0;
		if (cosPhi <= 0) {
			fr = SPECULAR_REDUCT * l->getRed() * surfaceMaterial.getSpecularRed() * phong;
			fg = SPECULAR_REDUCT * l->getGreen() * surfaceMaterial.getSpecularGreen() * phong;
			fb = SPECULAR_REDUCT * l->getBlue() * surfaceMaterial.getSpecularBlue() * phong;

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
		float cosTheta = inner_prod(surfaceNormal, l->getUnitVector()); 
		if (cosTheta >= 0) {
			fr =  DIFFUSE_FACT * surfaceMaterial.getDiffuseRed() * std::abs(cosTheta);
			fg =  DIFFUSE_FACT * surfaceMaterial.getDiffuseGreen() * std::abs(cosTheta);
			fb =  DIFFUSE_FACT * surfaceMaterial.getDiffuseBlue() * std::abs(cosTheta);
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
	int ar = surfaceMaterial.getAmbientRed() * AMB_LIGHT;
	int ag = surfaceMaterial.getAmbientGreen() * AMB_LIGHT;
	int ab = surfaceMaterial.getAmbientBlue() * AMB_LIGHT;

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
	ublas::vector<float> fc = ray.paraPos(far);	
	//fc (Z) = -c.getFarClip();
	float disN = norm_2(point - ray.getPRP());
	float disF = norm_2(fc - point);

	int pdepth = (int)  COLOR_MAX - std::min( float (COLOR_MAX), (COLOR_MAX) * (disN/(far - near)));

	ublas::vector<int> returnVal(VECTOR_C);
	returnVal(RED) = pr;
	returnVal(GREEN) = pg;
	returnVal(BLUE) = pb;
	returnVal(ALPHA) = pdepth;
	
	return returnVal;
}
