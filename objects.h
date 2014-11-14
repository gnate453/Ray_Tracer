#ifndef OBJS_H
#define OBJS_H 1

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <cmath>
#include <list>
#include <map>
#include <queue>

#ifndef NAME_UBLAS_ALIAS
#define NAME_UBLAS_ALIAS 1
namespace ublas = boost::numeric::ublas;
#endif /*namespace alias for boost uBLAS*/

#define VECTOR_2D 2
#define VECTOR_3D 3
#define VECTOR_2DH 3
#define VECTOR_3DH 4
#define	VECTOR_C 4
#define MATRIX_2D 2
#define MATRIX_3D 3
#define MATRIX_2DH 3
#define MATRIX_3DH 4
#define VERTEX_PER_FACE 3
#define P_ONE 0
#define P_TWO 1
#define P_THREE 2
#define X 0
#define Y 1
#define Z 2
#define W 3
#define RED 0
#define GREEN 1
#define BLUE 2
#define ALPHA 3
#define COLOR_MAX 255
#define COLOR_VALUES 256
#define AMB_LIGHT 20.0
#define DIFFUSE_FACT 0.00001
#define ZERO 0
#define BEGIN 0
#define LENGTH_1 1
#define LENGTH_2 2
#define LENGTH_6 6
#define LENGTH_7 7


ublas::vector<float> subtractVectors(ublas::vector<float> v1, ublas::vector<float> v2);

class Material {
	private:
	std::string name;
	ublas::vector<float> ka;
	ublas::vector<float> kd;
	ublas::vector<float> ks;

	public:
	Material();
	Material(const std::string&, const ublas::vector<float>&);
	Material(const std::string&, const ublas::vector<float>&, const ublas::vector<float>&, const ublas::vector<float>&);
	std::string getName();
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
	Face(const ublas::vector<float>&, const ublas::vector<float>&, const ublas::vector<float>& );
	ublas::vector<float> getVertex(int);
	ublas::vector<float> getNormal();
	bool isOnFace(ublas::vector<float>);
	
};

class Sphere {
	private:
	std::string name;
	Material color;	
	float radius;
	ublas::vector<float> originWorldCoord;

	public:
	Sphere(const std::string&, const Material&, const ublas::vector<float>&, const float&); 
	std::string getName();
	Material getColor(); 
	ublas::vector<float> getOrigin();
	float getRadius();
	float getRadiusSquared();
	float getDistanceToVPN(ublas::vector<float>);

};

class Polygon {
	private:
	std::string name;
	Material color;	
	std::list<Face> faces;

	public:
	Polygon(const std::string&, const Material&, std::list<Face>);
	std::string getName();
	Material getColor(); 
	std::list<Face> getFaces();
};

class Light {
	private:
	ublas::vector<float> directionVector;
	ublas::vector<int> color;

	public:
	Light(ublas::vector<float>, ublas::vector<float>);
	ublas::vector<float> getDirectionVector(); // v
	ublas::vector<float> getUnitVector();  // U = (v/||v||)
	float norm();	//||v|| = sqrt( (v1)^2 + (v2)^2 +...+(vn)^2 );
	ublas::vector<int> getColor();
	int getRed();
	int getGreen();
	int getBlue();
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
	std::list<Sphere> spheres;
	std::list<Polygon> polygons;
	std::list<Light> lights;
	std::map<std::size_t, ublas::vector<float> > vertices;
	std::map<std::string, Material> materials;

	public:
	std::list<Camera> getCameras();
	void addCamera(Camera);
	std::list<Scene> getScenes();
	void addScene(Scene);
	std::list<Sphere> getSpheres();
	void addSphere(Sphere);
	std::list<Polygon> getPolygons();
	void addPolygon(Polygon);
	std::list<Light> getLights();
	void addLight(Light);
	std::map<size_t, ublas::vector<float> > getVertices();
	void addVertex(ublas::vector<float>);
	std::map<std::string, Material> getMaterials();
	void addMaterial(std::string, Material);
};

#endif  /*define OBJS_H end */
