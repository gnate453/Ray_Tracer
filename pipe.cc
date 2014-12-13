#include "pipe.h"

float smoothDecimal(float d, int n) {
	return (std::floor(d * pow(BASE_TEN, n) + 0.5) / pow(BASE_TEN, n));
}

ublas::vector<float> smoothColor(ublas::vector<float> c) {
	c (RED) = std::min(float (COLOR_MAX), c (RED));
	c (GREEN) = std::min(float (COLOR_MAX), c (GREEN));
	c (BLUE) = std::min(float (COLOR_MAX), c (BLUE));

	return c;
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
					float beta = (2.0/(sit->getHeight()-1))*i - 1;
					float alpha = (2.0/(sit->getWidth()-1))*j - 1;
					
					ublas::vector<float> r(VECTOR_3D);
					ublas::vector<float> pixelWorldCoord (VECTOR_3D);

					r (X) = (-cit->getNearClip() * cit->getVPN()(X)) + (alpha * camU(X)) + (beta * camV(X));
					r (Y) = (-cit->getNearClip() * cit->getVPN()(Y)) + (alpha * camU(Y)) + (beta * camV(Y));
					r (Z) = (-cit->getNearClip() * cit->getVPN()(Z)) + (alpha * camU(Z)) + (beta * camV(Z));
	
					Ray tmpR(cit->getPRP(), r);
		
					//std::cout<<"\nPixel: "<<j<<", "<<i<<std::endl;
					ublas::vector<float> color = rayTrace(tmpR, &w, sit->getRecursionDepth(), true); 

					color = smoothColor(color);

					//std::cout<<"Final color: "<<color<<std::endl;

					img->setPixelRed(j, i, (int) color(RED));
					img->setPixelGreen(j, i, (int) color(GREEN));
					img->setPixelBlue(j, i,(int) color(BLUE));
					img->setPixelDepth(j, i,(int) color(ALPHA));


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

ublas::vector<float> rayTrace(Ray r, World *w, int depthCounter, bool f) {

	std::list<Light> lights = w->getLights();
	ublas::vector<float> tracedColor (VECTOR_C);
	tracedColor (RED) = 0.0;
	tracedColor (GREEN) = 0.0;
	tracedColor (BLUE) = 0.0;
	tracedColor (ALPHA) = 0.0;
	float cosThetaRay = 0;
	bool first = f;
		
	if (depthCounter >= 0) {
		//TODO intersect RAY GIVEN (r) with spheres and/or polygons.
		Intersection* closest = NULL;
		if (!w->getSpheres().empty() && !w->getPolygons().empty()) {
			Intersection* closestSphere = intersectRayWithSpheres(r, w->getSpheres());
			Intersection* closestPolygon = intersectRayWithPolygons(r, w->getPolygons());
			if (closestSphere != NULL && closestPolygon != NULL) {
				if (closestSphere->getDepth() < closestPolygon->getDepth()) {
					closest = closestSphere;
					delete closestPolygon;
				}
				else {
					closest = closestPolygon;
					delete closestSphere;
				}
			}
			else if (closestSphere != NULL && closestPolygon == NULL) {
				closest = closestSphere;
			}
			else if (closestSphere == NULL && closestPolygon != NULL) {
				closest = closestPolygon;
			}
		}
		else if (!w->getSpheres().empty() &&  w->getPolygons().empty()) {
			Intersection* closestSphere = intersectRayWithSpheres(r, w->getSpheres());
			closest = closestSphere;
		}
		else if (w->getSpheres().empty() &&  !w->getPolygons().empty()) {
			Intersection* closestPolygon = intersectRayWithPolygons(r, w->getPolygons());
			closest = closestPolygon;
		}
		else
			return tracedColor;
		
		
		//std::cout<<"Intersection: "<<closest->isReal()<<std::endl;
		
		//the given ray intersects with an object.
		if (closest != NULL) {
			
			//std::cout<<"r from: "<<r.getPRP()<<" U: "<<r.unitVector()<<std::endl;

			//calculate reflection of the given ray on closest surface
			cosThetaRay = inner_prod(r.unitVector(), closest->getSurfaceNormal());
			ublas::vector<float> R = (2 * cosThetaRay * closest->getSurfaceNormal()) - r.unitVector();			
			Ray newRay = Ray(closest->getPoint(), R);
			
			//std::cout<<"new orgin: "<<closest->getPoint()<<" R: "<<R<<std::endl;
			//std::cout<<"surface Normal: "<<closest->getSurfaceNormal()<<std::endl;

			//std::cout<<"distance: "<<closest->getDepth()<<std::endl;
			//std::cout<<"diffuse color: "<<closest->getSurfaceMaterial().getDiffuseProperties()<<std::endl;
			//std::cout<<"Specular color: "<<closest->getSurfaceMaterial().getSpecularProperties()<<std::endl;
			
			//recurse to the next ray.  Then recieve the color of it's specular highlighs
			tracedColor = rayTrace(newRay, w, depthCounter-1, false);
			
			//add your specular high lights.  
			for (std::list<Light>::iterator	l = lights.begin(); l != lights.end(); ++l) {
				bool blocked = false;	
							
				Ray L(closest->getPoint(), l->getDirectionVector(closest->getPoint()));
				//std::cout<<"for light in direction: "<<L.rayVector()<<std::endl;
				
				if (!w->getSpheres().empty()) {
 					Intersection* blocker = intersectRayWithSpheres(L, w->getSpheres());
					if (blocker != NULL) {
						blocked = true;
						//std::cout<<"blocked"<<std::endl;
						//std::cout<<"this point: "<<closest->getPoint()<<std::endl;
						//std::cout<<"blocking point"<<blocker->getPoint()<<std::endl;
						delete blocker;
					}
				}
				if (!blocked && !w->getPolygons().empty()) {
					Intersection* blocker = intersectRayWithPolygons(L, w->getPolygons());
					if (blocker != NULL) {
						blocked = true;
						delete blocker;
					}
				}

				//if (!blocked) {
					//calculate reflection of L
					 
					//std::cout<<"not blocked!!"<<std::endl;
					float cosThetaDiffuse = inner_prod(closest->getSurfaceNormal(), L.unitVector());
					float cosThetaLight = inner_prod(L.unitVector(), closest->getSurfaceNormal());
					ublas::vector<float> H = (2 * cosThetaLight * closest->getSurfaceNormal()) - L.unitVector() ;
					float cosPhiLight = inner_prod(r.unitVector(), H);
						
					//std::cout<<"reflection is H: "<<H<<std::endl;

					if (cosPhiLight <= 0) {
						tracedColor += calcSpecularColor(closest->getSurfaceMaterial(), l->getColor(), cosPhiLight);
						
						//std::cout<<"Added specular color"<<std::endl;
						//std::cout<<"color: "<<tracedColor<<std::endl;
					}
					if (first) {
						if (cosThetaDiffuse > 0) { 
							tracedColor += calcDiffuseColor(closest->getSurfaceMaterial(), l->getColor(), cosThetaDiffuse);
							//std::cout<<"Added Diffuse color, after returned to first hit!"<<std::endl;
							//std::cout<<"color: "<<tracedColor<<std::endl;
						}
					}
				//}
			}
			if (first) {
				tracedColor += calcAmbientColor(closest->getSurfaceMaterial());
			}

			delete closest;	 
		}
		else {
			//no intersection so return no color.
			return tracedColor;
		}
	}
	return tracedColor;
}

Intersection* intersectRayWithSpheres(Ray ray, std::list<Sphere> spheres) {

	ublas::vector<float> U = ray.unitVector();
	float closest = FAR_FAR_AWAY;
	Intersection* closestIntersection = NULL;
	//std::cout<<" U: "<<U<<std::endl;
	for (std::list<Sphere>::iterator s = spheres.begin(); s != spheres.end(); ++s)
	{
		float v, csqd, dsqd;
		
		v = inner_prod(s->getOrigin() - ray.getPRP(), U);
		csqd = s->getDistanceToPixel(ray.getPRP()) * s->getDistanceToPixel(ray.getPRP());
		dsqd = s->getRadiusSquared() - (csqd - (v * v));

		if (dsqd > 0) {
			float d = std::sqrt(dsqd);
			//std::cout<<"  d^2: "<<dsqd<<" d: "<<d<<"\n"<<std::endl;
			//normal to sphere
			ublas::vector<float> S = ray.paraPos(v-d);

			ublas::vector<float> N = S - s->getOrigin();
			N = (1/norm_2(N)) * N; 
				
			//record closest distance interestion with infromation to color pixel.
			if ((v-d) < closest) {
				closest = v-d;
				closestIntersection = new Intersection( (v-d), S, N, s->getColor(), true);
				//std::cout<<"In intersect distance: "<<closestIntersection->getDepth()<<std::endl;
				//std::cout<<"diffuse color: "<<closestIntersection->getSurfaceMaterial().getDiffuseProperties()<<std::endl;
				//std::cout<<"Specular color: "<<closestIntersection->getSurfaceMaterial().getSpecularProperties()<<std::endl;
			}	
		}	//end if intersection possible
	}	//end for each sphere

	return closestIntersection;
} 

Intersection* intersectRayWithPolygons(Ray ray, std::list<Polygon> polygons) {
	ublas::vector<float> rW = -ray.unitVector();
	ublas::vector<float> e1;
	ublas::vector<float> e2;
	ublas::vector<float> A;

	float closest = FAR_FAR_AWAY;
	Intersection* closestIntersection = NULL;

	for	(std::list<Polygon>::iterator p = polygons.begin(); p != polygons.end(); ++p) {
		float disToClosestFace = FAR_FAR_AWAY;
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
			float tstar = answer (Z);
				
			if (tstar < closest){
				closest = tstar;
				closestIntersection = new Intersection(tstar, ray.paraPos(tstar), f->getNormal(), p->getColor(), true);
			}
		} //end for each face
	} //end for each polygon

	return closestIntersection;
}

ublas::vector<float> calcSpecularColor(Material k, ublas::vector<float> b, float cosPhi) {

	ublas::vector<float> color (VECTOR_C);	
	float phong = pow(cosPhi, k.getSpecularAlpha());
	//phong = smoothDecimal(phong, COLOR_PERCISION);
	phong = std::abs(phong);
	//if (phong < EPSILON && phong > 0)
	//	phong = EPSILON;

	//std::cout<<"\n Specular calculator: "<<std::endl;
	//std::cout<<"light color: "<<b<<std::endl;
	//std::cout<<"Material Spec: "<<k.getSpecularProperties()<<std::endl;
	//std::cout<<"cosPhi: "<<cosPhi<<std::endl;
	//std::cout<<"phong: "<<phong<<std::endl;

	color (RED) = (b (RED) * k.getSpecularRed() * phong);
	color (GREEN) = (b (GREEN) * k.getSpecularGreen() * phong);
	color (BLUE) = (b (BLUE) * k.getSpecularBlue() * phong);

	//TODO: Scale by COLOR_MAX?
	//TODO: ALPHA?
	color(ALPHA) = 0.0;

	//std::cout<<"Specular calculator returning: "<<color<<std::endl;

	return color;
}

ublas::vector<float> calcDiffuseColor(Material k, ublas::vector<float> b, float cosTheta) {

	ublas::vector<float> color (VECTOR_C);
	
	//cosTheta = smoothDecimal(cosTheta, COLOR_PERCISION);
	
	//std::cout<<"\n Diffuse calculator: "<<std::endl;
	//std::cout<<"light color: "<<b<<std::endl;
	//std::cout<<"Material Diff: "<<k.getDiffuseProperties()<<std::endl;
	//std::cout<<"cosTheta: "<<cosTheta<<std::endl;

	//calc collor
	color (RED) = (b (RED) * k.getDiffuseRed() * cosTheta);
	color (GREEN) = (b (GREEN) * k.getDiffuseGreen() * cosTheta);
	color (BLUE) = (b (BLUE) * k.getDiffuseBlue() * cosTheta);

	//TODO: scale by COLOR_MAX?
	//TODO: ALPHA?
	color(ALPHA) = 0.0;

	
	//std::cout<<"Diffuse calculator returning: "<<color<<std::endl;

	return color;
}

ublas::vector<float> calcAmbientColor(Material k) {
	
	ublas::vector<float> color (VECTOR_C);

	color (RED) = k.getAmbientRed() * AMB_LIGHT; 
	color (GREEN) = k.getAmbientGreen() * AMB_LIGHT;
	color (BLUE) = k.getAmbientBlue() * AMB_LIGHT;

	//TODO: scale by COLOR_MAX?
	//TODO: ALPHA?
	color(ALPHA) = 0.0;

	return color;
}
