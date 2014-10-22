#include "objects.h"

Sphere::Sphere(std::string n, ublas::vector<float> o,
				ublas::vector<float> c, float r){
	name = n;
	originWorldCoord = o;
	color = c;
	radius = r;
} 

void Sphere::setRadiusSquared(float i) {
	radiusSquared = i;
}

void Sphere::setDistanceToPRPSquared(float i) {
	disToPRPSquared = i;
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

float Sphere::getRadiusSquared() {
	return radiusSquared;
}

float Sphere::getDistanceToPRPSquared() {
	return disToPRPSquared;
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


Ray::Ray(ublas::vector<float> prp) {
	focusWorldCoord = prp;
}

void Ray::setX(int x) {
	screenX = x;
}

int Ray::getX() {
	return screenX;
}

void Ray::setY(int y) {
	screenY = y;
}

int Ray::getY() {
	return screenY;
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
	
Image::Image(std::string n, int w, int h) 
{
	name = n;
	width = w;
	height = h;
	initImageData();
}

Image::Image(const Image& c) {
	name = c.getName();
	width = c.getWidth();
	height = c.getHeight();
	imgColorRed = new float [width*height];
	*imgColorRed = c.getImgRedData();
	imgColorGreen = new float [width*height];
	*imgColorGreen = c.getImgGreenData();
	imgColorBlue = new float [width*height];
	*imgColorBlue = c.getImgBlueData();
	imgDepthData = new float [width*height];
	*imgDepthData = c.getImgDepthData();
}

Image Image::operator=(Image rhs) {
	name = rhs.getName();
	width = rhs.getWidth();
	height = rhs.getHeight();
	imgColorRed = new float [width*height];
	*imgColorRed = rhs.getImgRedData();
	imgColorGreen = new float [width*height];
	*imgColorGreen = rhs.getImgGreenData();
	imgColorBlue = new float [width*height];
	*imgColorBlue = rhs.getImgBlueData();
	imgDepthData = new float[width*height];
	*imgDepthData = rhs.getImgDepthData();

	return *this; 
}

void Image::initImageData() {
	imgColorRed = new float [width*height];
	imgColorGreen = new float [width*height];
	imgColorBlue = new float [width*height];
	imgDepthData = new float [width*height];	
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

float Image::getImgRedData() const {
	return *imgColorRed;
}

float Image::getImgGreenData() const {
	return *imgColorGreen;
}

float Image::getImgBlueData() const {
	return *imgColorBlue;
}

float Image::getImgDepthData() const {
	return *imgDepthData;
}

void Image::setPixelRed(int x, int y, float p) {	
	imgColorRed[ (y*width) + x ] = p;
}

void Image::setPixelGreen(int x, int y, float p) {	
	imgColorGreen[ (y*width) + x ] = p;
}

void Image::setPixelBlue(int x, int y, float p) {	
	imgColorBlue[ (y*width) + x ] = p;
}

float Image::getPixelRed(int x,int y) {
	return imgColorRed[ (y*width)+ x ];
}

float Image::getPixelGreen(int x,int y) {
	return imgColorGreen[ (y*width)+ x ];
}

float Image::getPixelBlue(int x,int y) {
	return imgColorBlue[ (y*width)+ x ];
}


void Image::setPixelDepth(int x, int y, float d) {
	imgDepthData[ (y*width) + x ] = d;
}

float Image::getPixelDepth(int x,int y) {
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
