#include "objects.h"


Material::Material(const &ublas::vector<float> c) {
	ka = c;
	kd = c;
	ks = c;
}

Material::Material(const ublas::vector<float>& a, const ublas::vector<float>& d, const ublas::vector<float>& s) {
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
	reutn ka(GREEN);
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
	
ublas::vector<float> Face::getNormal() {
	e1 = p2 - p1;
	e2 = p3 - p2;
	
	return outer_ prod(e2, e1);
}

bool Face::isOnFace(ublas::vector<float> p) {
	if (inner_prod(getNormal(), p) == 0)
		return true;
	else 
		return false;	
}

Object::Object(const std::string& n, const Material& mtl) {
	
}

std::string Object::getName(){
	return name;
}

ublas::vector<float> Object::getColor(){
	return color;
}

void Object::setDistanceToPRPSquared(float i) {
	disToPRPSquared = i;
}

float Object::getDistanceToPRPSquared() {
	return disToPRPSquared;
}

Sphere::Sphere(const std::string& n, const Material& mtl, const ublas::vector<float>& o,const float &r) : Object(n, mtl) {
	originWorldCoord = o;
	radius = r;
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

Polygon::Polygon(const std::string& n, const Material& mtl, std::list<faces> f) : Object(n, mtl) {
	faces = f;
}

std::list<Face> getFaces() {
	return faces;
}

int getFaceCount() {
	return faces.size();
}

Light::Light( ublas::vector<float> dir, ublas vector<float> l) {
	directionVector = dir;
	color = l; 
}

ublas::vector<float> Light::getDirectionVector() {
	return directionVector;
}

ublas::vector<float> Light::getUnitVector() {
	return (1/norm()) * directionVector();
}

float Light::norm() {
	return norm_2(directionVector);
}

ublas::vector<float> Light::getColor() {
	return color();
}

float Light::getRed() {
	color(RED);
}

float Light::getGreen() {
	color(GREEN);
}

float Light::getBlue() {
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
	return getPixel() - getPRP(); 
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
	return ( getPixel() + unitVectorScaled(s) );
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

std::list<Object> World::getObjects() {
	return objects;
}

void World::addObject(Object o) {
	objects.push_back(o);
}

std::map<ublas::vector<float>> World::getVetices() {
	return verticies;
}

void World::addVertex(ublas::vector<float> v) {
	vertices [vertices.size()] = v;

std::map<ublas::vector<float>> World::getMaterials() {
	return materials;
}

void World::addMaterial(Material mtl) {
	materials [materials.size()] =  mtl;
