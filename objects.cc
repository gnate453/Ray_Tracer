#include "objects.h"

ublas::vector<float> crossProductVectors(ublas::vector<float> v1, ublas::vector<float> v2) {
	ublas::vector<float> r (VECTOR_3D);

	r (X) = ((v1(Y)) * v2(Z)) - ((v1(Z)) * v2(Y));
	r (Y) = ((v1(X)) * v2(Z)) - ((v1(Z)) * v2(X));
	r (Z) = ((v1(X)) * v2(Y)) - ((v1(Y)) * v2(X));

	return r;
}

Material::Material() {
	name = "Default";
	
	ublas::vector<float> a (VECTOR_C);
	a (RED) = 1.0;
	a (GREEN) = 1.0;
	a (BLUE) = 0.0;
	ka = a;	

	ublas::vector<float> d (VECTOR_C);
	d (RED) = 1.0;
	d (GREEN) = 0.0;
	d (BLUE) = 0.0;
	kd = d;	

	ublas::vector<float> s (VECTOR_C);
	s (RED) = 1.0;
	s (GREEN) = 1.0;
	s (BLUE) = 0.0;
	s (ALPHA) = 1.0;
	ks = s;	
}

Material::Material(const std::string& n, const ublas::vector<float>& c) {
	name = n;
	ka = c;
	kd = c;
	ks = c;
}

Material::Material(const std::string& n, const ublas::vector<float>& a, const ublas::vector<float>& d, const ublas::vector<float>& s) {
	name = n;
	ka = a;
	kd = d;
	ks = s;
}

ublas::vector<float> Material::getAmbientProperties() {
	return ka;
}

float Material::getAmbientRed() {
	return ka(RED);
}

float Material::getAmbientGreen() {
	return ka(GREEN);
}

float Material::getAmbientBlue() {
	return ka(BLUE);
}

ublas::vector<float> Material::getDiffuseProperties() {
	return kd;
}

float Material::getDiffuseRed() {
	return kd(RED);
}

float Material::getDiffuseGreen() {
	return kd(GREEN);
}

float Material::getDiffuseBlue() {
	return kd(BLUE);
}

ublas::vector<float> Material::getSpecularProperties() {
	return ks;
}

float Material::getSpecularRed() {
	return ks(RED);
}

float Material::getSpecularGreen() {
	return ks(GREEN);
}

float Material::getSpecularBlue() {
	return ks(BLUE);
}

float Material::getSpecularAlpha() {
	return ks(ALPHA);
}

Face::Face(const ublas::vector<float>& v1, 
			const ublas::vector<float>& v2, 
			const ublas::vector<float>& v3) {
	p1 = v1;
	p2 = v2;
	p3 = v3;
}

ublas::vector<float> Face::getVertex(int i) {
	switch (i) {
	case P_ONE :
		return p1;
	case P_TWO :
		return p2;
	case P_THREE :
		return p3;
	}
}
	
ublas::vector<float> Face::getNormal() {
	ublas::vector<float> e1 = p2 - p1;
	ublas::vector<float> e2 = p3 - p1;

	return (1/norm_2(crossProductVectors(e1, e2))) * crossProductVectors(e1, e2);
}

bool Face::isOnFace(ublas::vector<float> p) {
	if (inner_prod(getNormal(), p) == 0)
		return true;
	else 
		return false;	
}

Sphere::Sphere(const std::string& n, const Material& mtl, const ublas::vector<float>& o,const float &r) {
	name = n;
	color = mtl;
	originWorldCoord = o;
	radius = r;
} 

std::string Sphere::getName() {
	return name;
}

Material Sphere::getColor() {
	return color;
}

ublas::vector<float> Sphere::getOrigin(){
	return originWorldCoord;
}

float Sphere::getRadius() {
	return radius;
}

float Sphere::getRadiusSquared() {
	return radius * radius;
}

float Sphere::getDistanceToPixel(ublas::vector<float> pixel) {
	return norm_2(originWorldCoord - pixel);
}

//Polygon::Polygon(const Polygon& c) {
//	name = c.getName();
//	color = c.getColor();
//	faces = c.getFaces();
//}

Polygon::Polygon(const std::string& n, const Material& mtl, std::list<Face> f) {
	name = n;
	color = mtl;
	faces = f;
}

//Polygon Polygon::operator=(Polygon rhs) {
//	name = rhs.getName();
//	color = rhs.getColor();
//	faces = rhs.getFaces();
	
//	return *this;
//}

std::string Polygon::getName() {
	return name;
}

Material Polygon::getColor() {
	return color;
}

std::list<Face> Polygon::getFaces() {
	return faces;
}

Light::Light( ublas::vector<float> dir, ublas::vector<float> l) {
	directionVector = dir;
	color = l; 
}

ublas::vector<float> Light::getDirectionVector() {
	return directionVector;
}

ublas::vector<float> Light::getUnitVector() {
	return (1/norm()) * getDirectionVector();
}

float Light::norm() {
	return norm_2(getDirectionVector());
}

ublas::vector<int> Light::getColor() {
	return color;
}

int Light::getRed() {
	color(RED);
}

int Light::getGreen() {
	color(GREEN);
}

int Light::getBlue() {
	color(BLUE);
}

Ray::Ray(ublas::vector<float> prp, ublas::vector<float> rayVector, ublas::vector<float> vppc, int x, int y, float a, float b) {
	focusWorldCoord = prp;
	r = rayVector;
	pixelWorldCoord = vppc;
	screenX = x;
	screenY = y;
	alpha = a;
	beta = b;
}

int Ray::getScreenX() {
	return screenX;
}

int Ray::getScreenY() {
	return screenY;
}

float Ray::getAlpha() {
	return alpha;
}

float Ray::getBeta() {
	return beta;
}

ublas::vector<float> Ray::getPixelWorldCoord() {
	return pixelWorldCoord;
}

ublas::vector<float> Ray::getPRP() {
	return focusWorldCoord;			
}

//v = L - E, L is pixel of view plane, E is PRP
ublas::vector<float> Ray::rayVector() {
	return r; 
} 	

// ||v|| = sqrt( (v1)^2 + (v2)^2 +...+(vn)^2 )
float Ray::norm() {
	return norm_2(rayVector());
}

// U = (v/||v||)
ublas::vector<float> Ray::unitVector() {
	return ( (1/norm()) * rayVector() );
}		

//sU	
ublas::vector<float> Ray::unitVectorScaled(float s) {
	return (s * unitVector());
}	
	
// R(s) = L + sU,  L is pixel of view plane.
ublas::vector<float> Ray::paraPos(float s) {
	return getPRP() + unitVectorScaled(s);
}

Camera::Camera(std::string n, ublas::vector<float> prp,
				 ublas::vector<float> vpn, ublas::vector<float> vup,
				 float nc, float fc) {
	name = n;
	focusWorldCoord = prp;
	viewPlaneNormal = vpn;
	vectorUp = vup;
	nearClip = nc;
	farClip = fc;
	viewPlaneWorldCoord = prp - (nearClip * vpn);

	ublas::vector<float> nVPN = ((1/norm_2(vpn)) * vpn) ;
	horizontalAxis = (1/norm_2(crossProductVectors(vup, nVPN))) 
										* crossProductVectors(vup, nVPN);
	verticalAxis = crossProductVectors(nVPN, horizontalAxis);
}

std::string Camera::getName() {
	return name;
}
	
ublas::vector<float> Camera::getPRP() {
	return focusWorldCoord;
}

ublas::vector<float> Camera::getVRP() {
	return viewPlaneWorldCoord;
}

ublas::vector<float> Camera::getVPN() {
	return viewPlaneNormal;
}

ublas::vector<float> Camera::getVUP() {
	return vectorUp;
}

ublas::vector<float> Camera::getHorizontalVector() {
	return horizontalAxis;
}

ublas::vector<float> Camera::getVerticalVector() {
	return verticalAxis;
}

float Camera::getNearClip() {
	return nearClip;
}

float Camera::getFarClip() {
	return farClip;
}

Scene::Scene(std::string n, int w, int h, int depth) {
	name = n;
	width = w;
	height = h;
	recursionDepth = depth;
}

std::string Scene::getName() {
	return name;
}

int Scene::getWidth() {
	return width;
}

int Scene::getHeight(){
	return height;
}

int Scene::getRecursionDepth() {
	return recursionDepth;
}
	
Image::Image(std::string n, int w, int h) 
{
	name = n;
	width = w;
	height = h;
	initImageData();
	initPixelsWithData();
}

Image::Image(const Image& c) {
	name = c.getName();
	width = c.getWidth();
	height = c.getHeight();
	
	imgColorRed = new int [width*height];
	imgColorGreen = new int [width*height];
	imgColorBlue = new int [width*height];
	imgDepthData = new int [width*height];
	for (int i =0; i < (width*height); ++i)
	{
		imgColorRed[i] = c.getImgRedData()[i];
		imgColorGreen[i] = c.getImgGreenData()[i];
		imgColorBlue[i] = c.getImgBlueData()[i];
		imgDepthData[i] = c.getImgDepthData()[i];
	}

}

Image Image::operator=(Image rhs) {
	name = rhs.getName();
	width = rhs.getWidth();
	height = rhs.getHeight();
	
	imgColorRed = new int [width*height];
	imgColorGreen = new int [width*height];
	imgColorBlue = new int [width*height];
	imgDepthData = new int [width*height];
	
	for (int i =0; i < (width*height); ++i)
	{
		imgColorRed[i] = rhs.getImgRedData()[i];
		imgColorGreen[i] = rhs.getImgGreenData()[i];
		imgColorBlue[i] = rhs.getImgBlueData()[i];
		imgDepthData[i] = rhs.getImgDepthData()[i];
	}

	return *this; 
}

void Image::initImageData() {
	imgColorRed = new int [width*height];
	imgColorGreen = new int [width*height];
	imgColorBlue = new int [width*height];
	imgDepthData = new int [width*height];	
}

void Image::initPixelsWithData() {

	for (int i =0; i < (width*height); ++i)
	{
		imgColorRed[i] = ZERO;
		imgColorGreen[i] = ZERO;
		imgColorBlue[i] = ZERO;
		imgDepthData[i] = ZERO;
	}

}

std::string Image::getName() const {
	return name;
}

int Image::getWidth() const {
	return width;
}

int Image::getHeight() const {
	return height;
}

int* Image::getImgRedData() const {
	return imgColorRed;
}

int* Image::getImgGreenData() const {
	return imgColorGreen;
}

int* Image::getImgBlueData() const {
	return imgColorBlue;
}

int* Image::getImgDepthData() const {
	return imgDepthData;
}

void Image::setPixelRed(int x, int y, int p) {	
	imgColorRed[ (y*width) + x ] = p;
}

void Image::setPixelGreen(int x, int y, int p) {	
	imgColorGreen[ (y*width) + x ] = p;
}

void Image::setPixelBlue(int x, int y, int p) {	
	imgColorBlue[ (y*width) + x ] = p;
}

int Image::getPixelRed(int x,int y) {
	return imgColorRed[ (y*width)+ x ];
}

int Image::getPixelGreen(int x,int y) {
	return imgColorGreen[ (y*width)+ x ];
}

int Image::getPixelBlue(int x,int y) {
	return imgColorBlue[ (y*width)+ x ];
}


void Image::setPixelDepth(int x, int y, int d) {
	imgDepthData[ (y*width) + x ] = d;
}

int Image::getPixelDepth(int x,int y) {
	return imgDepthData[ (y*width) + x ];
}
void Image::cleanImage() {
	delete[] imgColorRed;
	delete[] imgColorGreen;
	delete[] imgColorBlue;
	delete[] imgDepthData;
}

void Image::resetImage() {
	cleanImage();
	initImageData();
}

Image::~Image() {
	cleanImage();
}

std::list<Camera> World::getCameras() {
	return cameras;
}

void World::addCamera(Camera c) {
	cameras.push_back(c);
}

std::list<Scene> World::getScenes() {
	return scenes;
}

void World::addScene(Scene s) {
	scenes.push_back(s);
}

std::list<Sphere> World::getSpheres() {
	return spheres;
}

void World::addSphere(Sphere s) {
	spheres.push_back(s);
}

std::list<Polygon> World::getPolygons() {
	return polygons;
}

void World::addPolygon(Polygon p) {
	polygons.push_back(p);
}

std::list<Light> World::getLights() {
	return lights;
}

void World::addLight(Light l) {
	lights.push_back(l);
}

std::map<size_t, ublas::vector<float> > World::getVertices() {
	return vertices;
}

void World::addVertex(ublas::vector<float> v) {
	vertices [vertices.size()] = v;
}

std::map<std::string, Material> World::getMaterials() {
	return materials;
}

void World::addMaterial(std::string n, Material mtl) {
	materials [n] =  mtl;
}

Intersection::Intersection(float d, ublas::vector<float> p, ublas::vector<float> n, Material m) {
	depth = d;
	point = p;
	normal = n;
	surfaceMaterial = m;
}

float Intersection::getDepth() {
	return depth;
}

ublas::vector<float> Intersection::getPoint() {
	return point;
}

ublas::vector<float> Intersection::getSurfaceNormal() {
	return normal;
}

Material Intersection::getSurfaceMaterial() {
	return surfaceMaterial;
}

