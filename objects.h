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
#define AMB_LIGHT 10.0
#define ZERO 0
#define TOSS_F 2
#define TOSS_S 6

class Material {
	private:
	ublas::vector<float> ka;
	ublas::vector<float> kd;
	ublas::vector<float> ks;

	public:
	Material(ublas::vector<float>);
	Material(ublas::vector<float>, ublas::vector<float>, ublas::vector<float>);
	ublas::vector<float> getAmbientProperties();
	float getAmbientRed();
	float getAmbientGreen();
	float getAmbientBlue();
	ublas::vector<float> getDiffuseProperties();
	float getDiffuseRed();
	float getDiffuseGreen();
	float getDiffuseBlue();
	ublas::vector<float> getSpecularProperties();
	float getSpecularRed();
	float getSpecularGreen();
	float getSpecularBlue();

};

class Face {
	private:
	ublas::vector<float> p1;
	ublas::vector<float> p2;
	ublas::vector<float> p3;

	public:
	Face(ublas::vector<float>, ublas::vector<float>, ublas::vector<float>);
	ublas::vector<float> getVertex(int);
	float getNormal();
	bool isOnFace(ublas::vector<float>);
	
};

class Object {
	private:	
	std::string name;
	Material color;	
	float disToPRPSquared;
	
	public:
	Object(std::string, Material);	
	void setDistanceToPRPSquared(float);
	std::string getName();
	Material getColor(); 
	float getDistanceToPRPSquared();
};

class Sphere : public Object {
	private:
	ublas::vector<float> originWorldCoord;
	float radius;	
	float radiusSquared;

	public:
	Sphere(std::string, Material, ublas::vector<float>, float); 
	void setRadiusSquared(float);
	ublas::vector<float> getOrigin();
	float getRadius();
	float getRadiusSquared();
};

class Polygon : public Object {
	std::list<Face> faces;

	public:
	Polygon(std::string, Material, std::list<faces>);
	std::list<Face> getFaces();
	Face getFace(int);
	float getNormaltoFace(int);
	int getFaceCount();
};

class Light {
	private:
	ublas::vector<float> directionVector;
	ublas::vector<float> color;

	public:
	Light(ublas::vector<float>, ublas::vector<float>)
	ublas::vector<float> getDirectionVector(); // v
	ublas::vector<float> getUnitVector();  // U = (v/||v||)
	float norm();	//||v|| = sqrt( (v1)^2 + (v2)^2 +...+(vn)^2 );
	ublas::vector<float> getColor();
	float getRed();
	float getGreen();
	float getBlue();
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
	void setScreenX(int);
	int getScreenX();
	void setScreenY(int);
	int getScreenY();
	//TODO: set world coord should take the vup, and use d, ax, by to calculte world coordinate.
	void setPixelWorldCoord(ublas::vector<float>);
	ublas::vector<float> getPixelWorldCoord();
	void setPRP(ublas::vector<float>);
	ublas::vector<float> getPRP();
	ublas::vector<float> rayVector(); //v = L - E, L is pixel of view plane, E is PRP
	float norm();// ||v|| = sqrt( (v1)^2 + (v2)^2 +...+(vn)^2 )
	ublas::vector<float> unitVector();		// U = (v/||v||)
	ublas::vector<float> unitVectorScaled(float);		// sU	
	ublas::vector<float> paraPos(float);// R(s) = L + sU,  L is pixel of view plane.

};

//TODO: functionality for cameras ofset form the origin.
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

class Image {
	private:
	std::string name;
	int width;
	int height;
	int* imgColorRed;
	int* imgColorGreen;
	int* imgColorBlue;
	int* imgDepthData;
	
	public:
	Image(std::string, int, int);
	Image(const Image& c);
	Image operator=(Image rhs);
	void initPixelsWithData();
	std::string getName() const;
	int getWidth() const;
	int getHeight() const;
	void initImageData();
	int* getImgRedData() const;
	int* getImgGreenData() const;
	int* getImgBlueData() const;
	int* getImgDepthData() const;
	void setPixelRed(int,int,int);
	int getPixelRed(int,int);
	void setPixelGreen(int,int,int);
	int getPixelGreen(int,int);
	void setPixelBlue(int,int,int);
	int getPixelBlue(int,int);
	void setPixelDepth(int,int,int);
	int getPixelDepth(int,int);
	void cleanImage();
	void resetImage();
	~Image();

};

class World {
	private:
	std::list<Camera> cameras;
	std::list<Scene> scenes;
	std::list<Object> objects;
	std::list<ublas::vector<float>> vertices;
	std::list<Material> materials;

	public:
	std::list<Camera> getCameras();
	void addCamera(Camera);
	std::list<Scene> getScenes();
	void addScene(Scene);
	std::list<Object> getObjects();
	void addObject(Object);
	ublas::vector<float> getVertex(int);
	std::list<ublas::vector<float>> getVertices();
	void addVertex(ublas::vector<float>);
	Material getMaterial(int);
	std::list<Material> getMaterials();
	void addMaterial(Material);
};

#endif  /*define OBJS_H end */
