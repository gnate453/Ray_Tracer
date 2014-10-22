#ifndef OBJS_H
#define OBJS_H 1

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <list>

#ifndef NAME_UBLAS_ALIAS
#define NAME_UBLAS_ALIAS 1
namespace ublas = boost::numeric::ublas;
#endif /*namespace alias for boost uBLAS*/

#define VECTOR_2D 2
#define VECTOR_3D 3
#define	VECTOR_C 4
#define MATRIX_2D 3
#define MATRIX_3D 4
#define X 0
#define Y 1
#define Z 2
#define RED 0
#define GREEN 1
#define BLUE 2
#define ALPHA 3
#define COLOR_MAX 255
#define COLOR_VALUES 256

class Sphere {
	private:
	std::string name;
	ublas::vector<float> originWorldCoord;
	ublas::vector<float> color;	
	float radius;	
	float radiusSquared;
	float disToPRPSquared;

	public:
	Sphere(std::string, ublas::vector<float>, ublas::vector<float>, float); 
	void setRadiusSquared(float);
	void setDistanceToPRPSquared(float);
	std::string getName();
	ublas::vector<float> getOrigin();
	ublas::vector<float> getColor(); 
	float getRadius();
	float getRadiusSquared();
	float getDistanceToPRPSquared();

};

class Camera {
	private:
	std::string name;
	ublas::vector<float> focusWorldCoord;
	ublas::vector<float> viewPlaneWorldCoord;
	float nearClip;
	float farClip;
	
	public:
	Camera(std::string, ublas::vector<float>, ublas::vector<float>, float, float);
	std::string getName();
	ublas::vector<float> getPRP();
	ublas::vector<float> getVPN();
	float getNearClip();
	float getFarClip();

};

class Scene {
	private:	
	std::string name;
	int width;
	int height;
	int recursionDepth;	

	public:
	Scene(std::string, int, int, int);
	std::string getName();
	int getWidth();
	int getHeight();
	int getRecursionDepth();

};

class Ray {
	private:
	//world Coordinate of the Pixel this ray colors.
	ublas::vector<float> pixelWorldCoord;
	//world Coordinate of the PRP or foucas point.		
	ublas::vector<float> focusWorldCoord;
	int screenX;
	int screenY;		

	public:
	Ray(ublas::vector<float>);
	
	void setX(int);
	int getX();
	void setY(int);
	int getY();
	void setPixel(ublas::vector<float>);
	ublas::vector<float> getPixel();
	void setPRP(ublas::vector<float>);
	ublas::vector<float> getPRP();
	ublas::vector<float> rayVector(); //v = L - E, L is pixel of view plane, E is PRP
	float norm();// ||v|| = sqrt( (v1)^2 + (v2)^2 +...+(vn)^2 )
	ublas::vector<float> unitVector();		// U = (v/||v||)
	ublas::vector<float> unitVectorScaled(float);		// sU	
	ublas::vector<float> paraPos(float);// R(s) = L + sU,  L is pixel of view plane.

};

class World {
	private:
	std::list<Sphere> spheres;
	std::list<Camera> cameras;
	std::list<Scene> scenes;

	public:
	std::list<Sphere> getSpheres();
	void addSphere(Sphere);
	std::list<Camera> getCameras();
	void addCamera(Camera);
	std::list<Scene> getScenes();
	void addScene(Scene);
};

class Image {
	private:
	std::string name;
	int width;
	int height;
	float* imgColorRed;
	float* imgColorGreen;
	float* imgColorBlue;
	float* imgDepthData;
	
	public:
	Image(std::string, int, int);
	Image(const Image& c);
	Image operator=(Image rhs);
	std::string getName() const;
	int getWidth() const;
	int getHeight() const;
	void initImageData();
	float getImgRedData() const;
	float getImgGreenData() const;
	float getImgBlueData() const;
	float getImgDepthData() const;
	void setPixelRed(int,int,float);
	float getPixelRed(int,int);
	void setPixelGreen(int,int,float);
	float getPixelGreen(int,int);
	void setPixelBlue(int,int,float);
	float getPixelBlue(int,int);
	void setPixelDepth(int,int,float);
	float getPixelDepth(int,int);
	void cleanImage();
	void resetImage();
	~Image();

};
#endif  /*define OBJS_H end */
