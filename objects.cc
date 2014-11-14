#include "objects.h"

ublas::vector<float> subtractVectors(ublas::vector<float> v1, ublas::vector<float> v2) {
	if (v1.size() == VECTOR_3D && v2.size() == VECTOR_3D) {
		return v1 - v2;
	}
	else if (v1.size() == VECTOR_3DH && v2.size() == VECTOR_3D) {
		ublas::vector<float> r (VECTOR_3DH);
		r (X) = (v1 (X) / v1 (W)) - v2 (X);
		r (Y) = (v1 (Y) / v1 (W)) - v2 (Y);
		r (Z) = (v1 (Z) / v1 (W)) - v2 (Z);
		r (W) = 1.0; 
		return r;
	}
	else if (v1.size() == VECTOR_3D && v2.size() == VECTOR_3DH) {
		ublas::vector<float> r (VECTOR_3DH);
		r (X) = v1 (X) - (v2 (X) / v2 (W));
		r (Y) = v1 (Y) - (v2 (Y) / v2 (W));
		r (Z) = v1 (Z) - (v2 (Z) / v2 (W));
		r (W) = 1.0; 
		return r;
	}
	else if (v1.size() == VECTOR_3DH && v2.size() == VECTOR_3DH) {
		ublas::vector<float> r (VECTOR_3DH);
		r (X) = (v1 (X) / v1 (W)) - (v2 (X) / v2 (W));
		r (Y) = (v1 (Y) / v1 (W)) - (v2 (Y) / v2 (W));
		r (Z) = (v1 (Z) / v1 (W)) - (v2 (Z) / v2 (W));
		r (W) = 1.0; 
		return r;
	}

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
	d (GREEN) = 1.0;
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
	ublas::vector<float> e2 = p3 - p2;

	ublas::vector<float> n (VECTOR_3DH);

	n (X) = ((e2(Y)/e2(W)) * (e1(Z)/e1(W)))	- ((e2(Z)/e2(W)) * (e1(Y)/e1(W)));
	n (Y) = ((e2(X)/e2(W)) * (e1(Z)/e1(W)))	- ((e2(Z)/e2(W)) * (e1(X)/e1(W)));
	n (Z) = ((e2(X)/e2(W)) * (e1(Y)/e1(W)))	- ((e2(Y)/e2(W)) * (e1(X)/e1(W)));
	n (W) = 1;

	return n;
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

float Sphere::getDistanceToVPN(ublas::vector<float> vpn) {
	ublas::vector<float> t;
	if (vpn (Z) > originWorldCoord (Z)) {
		t = subtractVectors(originWorldCoord, vpn);
	}
	else {
		t = subtractVectors(vpn, originWorldCoord);
	}
	return norm_2(t);
}

Polygon::Polygon(const std::string& n, const Material& mtl, std::list<Face> f) {
	name = n;
	color = mtl;
	faces = f;
}

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
	return (1/norm()) * directionVector;
}

float Light::norm() {
	return norm_2(directionVector);
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

Ray::Ray(ublas::vector<float> prp) {
	focusWorldCoord = prp;
}

void Ray::setScreenX(int x) {
	screenX = x;
}

int Ray::getScreenX() {
	return screenX;
}

void Ray::setScreenY(int y) {
	screenY = y;
}

int Ray::getScreenY() {
	return screenY;
}

//TODO: this should take different params see comment in header.
void Ray::setPixelWorldCoord(ublas::vector<float> l) {
	pixelWorldCoord = l;
}

ublas::vector<float> Ray::getPixelWorldCoord() {
	return pixelWorldCoord;
}

void Ray::setPRP(ublas::vector<float> prp) {
	focusWorldCoord = prp;
}

ublas::vector<float> Ray::getPRP() {
	return focusWorldCoord;			
}

//v = L - E, L is pixel of view plane, E is PRP
ublas::vector<float> Ray::rayVector() {
	return getPixelWorldCoord() - getPRP(); 
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
	return ( getPixelWorldCoord() + unitVectorScaled(s) );
}

Camera::Camera(std::string n, ublas::vector<float> prp,
				 ublas::vector<float> vpn, float nc, float fc) {
	name = n;
	focusWorldCoord = prp;
	viewPlaneWorldCoord = vpn;
	nearClip = nc;
	farClip = fc;
}

std::string Camera::getName() {
	return name;
}
	
ublas::vector<float> Camera::getPRP() {
	return focusWorldCoord;
}

ublas::vector<float> Camera::getVPN() {
	return viewPlaneWorldCoord;
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
	for (int i =0; i < (width*height); ++i)
	{
		imgColorRed[i] = c.getImgRedData()[i];
	}
	imgColorGreen = new int [width*height];

	for (int i =0; i < (width*height); ++i)
	{
		imgColorGreen[i] = c.getImgGreenData()[i];
	}
	
	imgColorBlue = new int [width*height];

	for (int i =0; i < (width*height); ++i)
	{
		imgColorBlue[i] = c.getImgBlueData()[i];
	}
	
	imgDepthData = new int [width*height];
	
	for (int i =0; i < (width*height); ++i)
	{
		imgDepthData[i] = c.getImgDepthData()[i];
	}
}

Image Image::operator=(Image rhs) {
	name = rhs.getName();
	width = rhs.getWidth();
	height = rhs.getHeight();
	imgColorRed = new int [width*height];
	for (int i =0; i < (width*height); ++i)
	{
		imgColorRed[i] = rhs.getImgRedData()[i];
	}
	imgColorGreen = new int [width*height];

	for (int i =0; i < (width*height); ++i)
	{
		imgColorGreen[i] = rhs.getImgGreenData()[i];
	}
	
	imgColorBlue = new int [width*height];

	for (int i =0; i < (width*height); ++i)
	{
		imgColorBlue[i] = rhs.getImgBlueData()[i];
	}
	
	imgDepthData = new int [width*height];
	
	for (int i =0; i < (width*height); ++i)
	{
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
