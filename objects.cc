#include "objects.h"

Sphere::Sphere(std::string n, ublas::vector<float> o,
				ublas::vector<float> c, float r){
	name = n;
	originWorldCoord = o;
	color = c;
	radius = r;
} 

std::string Sphere::getName(){
	return name;
}

ublas::vector<float> Sphere::getOrigin(){
	return originWorldCoord;
}

ublas::vector<float> Sphere::getColor(){
	return color;
}

float Sphere::getRadius() {
	return radius;
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


Ray::Ray(ublas::vector<float> prp, ublas::vector<float> pixel) {
	pixelWorldCoord = pixel;
	focusWorldCoord = prp;
}

void Ray::setPixel(ublas::vector<float> l) {
	pixelWorldCoord = l;
}

ublas::vector<float> Ray::getPixel() {
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
	return getPixel()-getPRP(); 
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

std::list<Sphere> World::getSpheres() {
	return spheres;
}

void World::addSphere(Sphere s) {
	spheres.push_back(s);
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
	
Image::Image(int w, int h) 
{
	width = w;
	height = h;
	initImageData();
}

void Image::initImageData() {
	imgColorData = new ublas::vector<float> [width*height];
	imgDepthData = new float [width*height];	
}

void Image::setPixelColor(int x, int y, ublas::vector<float> p) {	
	imgColorData[ (y*width) + x ] = p;
}

ublas::vector<float> Image::getPixelColor(int x,int y) {
	return imgColorData[ (y*width)+ x ];
}


void Image::setPixelDepth(int x, int y, float d) {
	imgDepthData[ (y*width) + x ] = d;
}

float Image::getPixelDepth(int x,int y) {
	return imgDepthData[ (y*width) + x ];
}
void Image::cleanImage() {
	delete[] imgColorData;
	delete[] imgDepthData;
}

void Image::resetImage() {
	cleanImage();
	initImageData();
}

Image::~Image() {
	cleanImage();
}
