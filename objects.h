#ifndef OBJS_H
#define OBJS_H 1

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <cmath>
#include <list>
#include <map>
#include <queue>
#include <limits.h>

#ifndef NAME_UBLAS_ALIAS
#define NAME_UBLAS_ALIAS 1
namespace ublas = boost::numeric::ublas;
#endif /*namespace alias for boost uBLAS*/

#define VECTOR_2D 2
#define VECTOR_3D 3
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
#define SPHERE true
#define POLYGON false
#define AMB_LIGHT 20.0
#define LIGHT_FACT 20
#define EPSILON 0.000001
#define ZERO 0
#define ZERO_F 0.0
#define BEGIN 0
#define LENGTH_1 1
#define LENGTH_2 2
#define LENGTH_6 6
#define LENGTH_7 7


ublas::vector<double> crossProductVectors(ublas::vector<double> v1, ublas::vector<double> v2);


class Material {
	private:
	std::string name;
	ublas::vector<double> ka;
	ublas::vector<double> kd;
	ublas::vector<double> ks;

	public:
	Material();
	Material(const std::string&, const ublas::vector<double>&);
	Material(const std::string&, const ublas::vector<double>&, const ublas::vector<double>&, const ublas::vector<double>&);
	std::string getName();
	ublas::vector<double> getAmbientProperties();
	double getAmbientRed();
	double getAmbientGreen();
	double getAmbientBlue();
	ublas::vector<double> getDiffuseProperties();
	double getDiffuseRed();
	double getDiffuseGreen();
	double getDiffuseBlue();
	ublas::vector<double> getSpecularProperties();
	double getSpecularRed();
	double getSpecularGreen();
	double getSpecularBlue();
	double getSpecularAlpha();
};

class Face {
	private:
	ublas::vector<double> p1;
	ublas::vector<double> p2;
	ublas::vector<double> p3;

	public:
	Face(const ublas::vector<double>&, const ublas::vector<double>&, const ublas::vector<double>& );
	ublas::vector<double> getVertex(int);
	ublas::vector<double> getNormal();
	bool isOnFace(ublas::vector<double>);
	
};

class Sphere {
	private:
	std::string name;
	Material color;	
	double radius;
	ublas::vector<double> originWorldCoord;

	public:
	Sphere(const std::string&, const Material&, const ublas::vector<double>&, const double&); 
	std::string getName();
	Material getColor(); 
	ublas::vector<double> getOrigin();
	double getRadius();
	double getRadiusSquared();
	double getDistanceToPixel(ublas::vector<double>);

};

class Polygon {
	private:
	std::string name;
	Material color;	
	std::list<Face> faces;

	public:
	//Polygon(const Polygon&);
	Polygon(const std::string&, const Material&, std::list<Face>);
	//Polygon operator=(Polygon);
	//void copyFaces(std::list<Face>);
	std::string getName();	
	Material getColor(); 
	std::list<Face> getFaces();
	
};

class Light {
	private:
	ublas::vector<double> directionVector;
	ublas::vector<int> color;

	public:
	Light(ublas::vector<double>, ublas::vector<double>);
	ublas::vector<double> getDirectionVector(ublas::vector<double>); // v
	ublas::vector<double> getUnitVector(ublas::vector<double>);  // U = (v/||v||)
	double norm(ublas::vector<double>);	//||v|| = sqrt( (v1)^2 + (v2)^2 +...+(vn)^2 );
	ublas::vector<int> getColor();
	int getRed();
	int getGreen();
	int getBlue();
};

class Ray {
	private:
	//world Coordinate of the Pixel this ray colors. (L)
	ublas::vector<double> pixelWorldCoord;
	//direction vector of the ray;
	ublas::vector<double> r;
	//world Coordinate of the PRP or foucas point. (E)	
	ublas::vector<double> focusWorldCoord;
	int screenX;
	int screenY;
	double alpha;
	double beta;

	public:
	Ray(ublas::vector<double>, ublas::vector<double>, ublas::vector<double>, int, int, double, double);
	int getScreenX();
	int getScreenY();
	double getAlpha();
	double getBeta();
	//TODO: set world coord should take the vup, and use d, ax, by to calculte world coordinate.
	ublas::vector<double> getPixelWorldCoord();
	ublas::vector<double> getPRP();
	ublas::vector<double> rayVector(); //v = L - E, L is pixel of view plane, E is PRP
	double norm();// ||v|| = sqrt( (v1)^2 + (v2)^2 +...+(vn)^2 )
	ublas::vector<double> unitVector();		// U = (v/||v||)
	ublas::vector<double> unitVectorScaled(double);		// sU	
	ublas::vector<double> paraPos(double);// R(s) = L + sU,  L is pixel of view plane.

};

//TODO: functionality for cameras ofset form the origin.
class Camera {
	private:
	std::string name;
	ublas::vector<double> focusWorldCoord;
	ublas::vector<double> viewPlaneWorldCoord;
	ublas::vector<double> vectorUp;
	ublas::vector<double> viewPlaneNormal;
	ublas::vector<double> horizontalAxis;
	ublas::vector<double> verticalAxis;
	double nearClip;
	double farClip;
	
	public:
	Camera(std::string, ublas::vector<double>, ublas::vector<double>, ublas::vector<double>, double, double);
	std::string getName();
	ublas::vector<double> getPRP();
	ublas::vector<double> getVRP();
	ublas::vector<double> getVPN();
	ublas::vector<double> getVUP();
	ublas::vector<double> getHorizontalVector();
	ublas::vector<double> getVerticalVector();
	double getNearClip();
	double getFarClip();

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
	std::map<std::size_t, ublas::vector<double> > vertices;
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
	std::map<size_t, ublas::vector<double> > getVertices();
	void addVertex(ublas::vector<double>);
	std::map<std::string, Material> getMaterials();
	void addMaterial(std::string, Material);
};

class Intersection {
	private:
	double depth;
	Material surfaceMaterial;
	ublas::vector<double> point;
	ublas::vector<double> normal;

	public:
	Intersection(double, ublas::vector<double>, ublas::vector<double>, Material);
	double getDepth();
	ublas::vector<double> getPoint();
	ublas::vector<double> getSurfaceNormal();
	Material getSurfaceMaterial();
};

#endif  /*define OBJS_H end */
