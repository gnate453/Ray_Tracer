#include "pipe.h"

float smoothDecimal(float d, int n) {
	return (std::floor(d * pow(BASE_TEN, n) + 0.5) / pow(BASE_TEN, n));
}

World worldFromString(std::string d) {
	
	World newWorld;
	std::list<Face> facesInGroup;
	std::string currentGroupName = "Default";
	std::string curMaterial;
	std::string polyMaterial;
	std::stringstream data(d);
	std::string tmp;
	bool facesAdded = false;
	bool polygonAdded = false;
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

			if (facesAdded){		
				Polygon p(currentGroupName, (newWorld.getMaterials())[curMaterial], facesInGroup);
				polygonAdded = true;
				newWorld.addPolygon(p);
				facesInGroup.clear();
			}

			parse>>currentGroupName;
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

			//round diffuse and specular vales to nice decimal value!
			kd(RED) = smoothDecimal(kd(RED), COLOR_PERCISION);
			kd(GREEN) = smoothDecimal(kd(GREEN), COLOR_PERCISION);
			kd(BLUE) = smoothDecimal(kd(BLUE), COLOR_PERCISION);
			
			ks(RED) = smoothDecimal(ks(RED), COLOR_PERCISION);
			ks(GREEN) = smoothDecimal(ks(GREEN), COLOR_PERCISION);
			ks(BLUE) = smoothDecimal(ks(BLUE), COLOR_PERCISION);
			
			Material mtl(name, ka, kd, ks);
			newWorld.addMaterial(name, mtl);
		}
		//use material
		else if (*(tmp.begin()) == 'u') {	
			std::stringstream parse(tmp.substr(LENGTH_6, tmp.length()));
			//std::cout<<"Before use: "<<curMaterial<<std::endl;
			if (!polygonAdded)
				polyMaterial = curMaterial;

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
			//TODO: float CHECK this is still all that is needed for a scene.
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

			//std::cout<<"light color: "<<c<<std::endl;						

			newWorld.addLight(l);
		}
	}

	if ( (facesAdded && !groupNamed) || (facesAdded && groupNamed) ) {
		if (!facesInGroup.empty()) {
			if (polygonAdded) {
				Polygon p(currentGroupName, (newWorld.getMaterials())[curMaterial], facesInGroup);
				newWorld.addPolygon(p);
			}
			else {
				Polygon p(currentGroupName, (newWorld.getMaterials())[polyMaterial], facesInGroup);
				newWorld.addPolygon(p);
			}
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

			std::cout<<"U: "<<camU<<" V: "<<camV<<" PRP: "<<cit->getPRP()<<std::endl;

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
		if (dsqd > 0) {
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
	//ublas::matrix<float> M (VECTOR_3D, VECTOR_3D);
	//ublas::matrix<float> Mi (VECTOR_3D, VECTOR_3D);
	//ublas::matrix<float> Mb (VECTOR_3D, VECTOR_3D);
	//ublas::matrix<float> Mg (VECTOR_3D, VECTOR_3D);
	//ublas::matrix<float> Mt (VECTOR_3D, VECTOR_3D);
	ublas::vector<float> rW = -ray.unitVector();
	ublas::vector<float> e1;
	ublas::vector<float> e2;
	ublas::vector<float> A;

	float closest = c.getFarClip();
	Intersection closestIntersection = Intersection(closest, ray.paraPos(closest), rW, polygons.front().getColor());

	for	(std::list<Polygon>::iterator p = polygons.begin(); p != polygons.end(); ++p) {
		float disToClosestFace = c.getFarClip();
		std::list<Face> facesInP = p->getFaces();
		for (std::list<Face>::iterator f = facesInP.begin(); f != facesInP.end(); ++f) {
			e1  = f->getVertex(P_TWO) - f->getVertex(P_ONE);
			e2 = f->getVertex(P_THREE) - f->getVertex(P_ONE);
			A = ray.getPRP() - f->getVertex(P_ONE);
			
			//compute the det of [(p2-p1), (p3-p1), A] matrix
			float det = (e1(X) * ( (e2(Y)*rW(Z)) - (rW(Y)*e2(Z)) ))
				- (e2(X) * ( (rW(Z)*e1(Y)) - (rW(Y)*e1(Z)) ))
				- (rW(X) * ( (e1(Y)*e2(Z)) - (e2(Y)*e1(Z)) ));  	
			
			//if det is close to zero stop
			if (det == 0.0)
				continue;
			
			float inv_det = 1.0 / det;					

			float Ai = inv_det * ((e2(Y)*rW(Z)) - (rW(Y)*e2(Z)));
			float Bi = inv_det * -((e1(Y)*rW(Z)) - (rW(Y)*e1(Z)));
			float Ci = inv_det * ((e1(Y)*e2(Z)) - (e2(Y)*e1(Z)));
			float Di = inv_det * -((e2(X)*rW(Z)) - (rW(X)*e2(Z)));
			float Ei = inv_det * ((e1(X)*rW(Z)) - (rW(X)*e1(Z)));
			float Fi = inv_det * -((e1(X)*e2(Z)) - (e2(X)*e1(Z)));
			float Gi = inv_det * ((e2(X)*rW(Y)) - (rW(X)*e2(Y)));
			float Hi = inv_det * -((e1(X)*rW(Y)) - (rW(X)*e1(Y)));
			float Ii = inv_det * ((e1(X)*e2(Y)) - (e2(X)*e1(Y)));

			//result is the beta gamma tstar vector
			//ublas::vector<float> result (VECTOR_3D);
			ublas::vector<float> answer (VECTOR_3D);

			answer (X) = (Ai * A(X)) + (Di * A(Y)) + (Gi * A(Z));

			if (answer (X) < 0.0 || answer (X) > 1.0)
				continue;

			answer (Y) = (Bi * A(X)) + (Ei * A(Y)) + (Hi * A(Z));

			if (answer (Y) < 0.0 || answer (X) + answer (Y) > 1.0)				
				continue;

			answer (Z) = (Ci * A(X)) + (Fi * A(Y)) + (Ii * A(Z));

			//std::cout<<"result: "<<result<<"answer: "<<answer<<std::endl;

			if ( answer (Z) > c.getNearClip() && answer (Z) < c.getFarClip()  ) { 
				
				float tstar = answer (Z);
				
				if (tstar < closest){
					closest = tstar;
					closestIntersection = Intersection(closest, ray.paraPos(tstar), f->getNormal(), p->getColor());
				}
			} //end if intersection
		} //end for each face
	} //end for each polygon

	return closestIntersection;
}

ublas::vector<float> calcPixelColor(Ray ray, ublas::vector<float> point, ublas::vector<float> surfaceNormal, Material surfaceMaterial, std::list<Light> lights, float near, float far) {

	float intensity;
	float sr = 0.0;
	float sg = 0.0;
	float sb = 0.0;
	float dr = 0.0;
	float dg = 0.0;
	float db = 0.0;
	unsigned int pr = 0;
	unsigned int pg = 0;
	unsigned int pb = 0;
	int illuminationCountD = 0;
	int illuminationCountS = 0;
	for (std::list<Light>::iterator	l = lights.begin(); l != lights.end(); ++l) {
		
		//std::cout<<"Light: "<<l->getColor()<<std::endl;

		//Specular lighting.
		float cosPhi = inner_prod(ray.unitVector(), ((2 * inner_prod(l->getUnitVector(point),
									 surfaceNormal) * surfaceNormal) - l->getUnitVector(point)));	

		if (cosPhi <= 0) {
			float phong = pow(cosPhi, surfaceMaterial.getSpecularAlpha());
			smoothDecimal(phong, COLOR_PERCISION);
			phong = std::abs(phong);
			sr = sr + (l->getRed() * surfaceMaterial.getSpecularRed() * phong) / COLOR_MAX;
			sg = sg + (l->getGreen() * surfaceMaterial.getSpecularGreen() * phong) / COLOR_MAX;
			sb = sb + (l->getBlue() * surfaceMaterial.getSpecularBlue() * phong) / COLOR_MAX;

			illuminationCountS += 1;
		}

		//Diffuse Lighting
		float cosTheta = inner_prod(surfaceNormal, l->getUnitVector(point)); 
						// (norm_2(surfaceNormal) * norm_2(l->getUnitVector(point))) ; 

		//std::cout<<"cosTheta: "<<cosTheta<<std::endl;
		//std::cout<<"Object: "<<surfaceMaterial.getDiffuseProperties()<<std::endl;

		if (cosTheta > 0) {

			smoothDecimal(cosTheta, COLOR_PERCISION);
			dr = dr + (l->getRed() * surfaceMaterial.getDiffuseRed() * cosTheta) / COLOR_MAX;
			dg = dg + (l->getGreen() * surfaceMaterial.getDiffuseGreen() * cosTheta) / COLOR_MAX;
			db = db + (l->getBlue() * surfaceMaterial.getDiffuseBlue() * cosTheta) / COLOR_MAX;
			illuminationCountD +=1;
		}
	}			
	//std::cout<<"dr: "<<dr<<" : "<<dr/(illuminationCountD*COLOR_MAX)<<std::endl;
	//std::cout<<"dg: "<<dg<<" : "<<dg/(illuminationCountD*COLOR_MAX)<<std::endl;
	//std::cout<<"db: "<<db<<" : "<<db/(illuminationCountD*COLOR_MAX)<<std::endl;

	//Ambient Lighting.
	float ar = (surfaceMaterial.getAmbientRed() * AMB_LIGHT ) / COLOR_MAX;
	float ag = (surfaceMaterial.getAmbientGreen() * AMB_LIGHT / COLOR_MAX);
	float ab = (surfaceMaterial.getAmbientBlue() * AMB_LIGHT / COLOR_MAX);
	
	intensity = (sr+dr+ar) * LIGHT_FACT;
	pr = std::min(COLOR_MAX, (int) std::floor(COLOR_MAX * intensity));
	intensity = (sg+dg+ag) * LIGHT_FACT;
	pg = std::min(COLOR_MAX, (int) std::floor(COLOR_MAX * intensity));
	intensity = (sb+db+ab) * LIGHT_FACT;
	pb = std::min(COLOR_MAX, (int) std::floor(COLOR_MAX * intensity));

	//std::cout<<"red: s:"<<sr<<" d:"<<dr<<" i:"<<intensity<<std::endl;
	//std::cout<<"green: s:"<<sg<<" d:"<<dg<<" i:"<<intensity<<std::endl;
	//std::cout<<"blue: s:"<<sb<<" d:"<<db<<" i:"<<intensity<<std::endl;

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
