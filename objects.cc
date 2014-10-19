#include "objects.h"

Sphere::Sphere(std::string n, ublas::vector<float> o, ublas::vector<float> c, float r){
} 

std::string Sphere::getName(){
}

ublas::vector<float> Sphere::getOrigin(){
}

ublas::vector<float> Sphere::getColor(){
}

float Sphere::getRadius() {
}

Camera::Camera(std::string n, ublas::vector<float> prp, ublas::vector<float> vpn, float nc, float fc) {}
std::string Camera::getName() {}
ublas::vector<float> Camera::getPRP() {}
ublas::vector<float> Camera::getVPN() {}
float Camera::getNearClip() {}
float Camera::getFarClip() {}

Scene::Scene(std::string name, int w, int h, int depth) {}
std::string Scene::getName() {}
int Scene::getWidth() {}
int Scene::getHeight(){}
int Scene::getRecursionDepth() {}


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
	
Image::Image(int, int) {}
void Image::setPixel(int,int,ublas::vector<float>) {}
ublas::vector<float> Image::getPixel(int,int) {}
void Image::clear() {}
void Image::clean() {}
