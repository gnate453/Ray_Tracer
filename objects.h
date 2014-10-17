#ifndef OBJS_H
#define OBJS_H 1

#ifndef INC_BOOST_VECTOR
#define INC_BOOST_VECTOR 1
#include <boost/numeric/ublas/vector.hpp>
#endif  /*include boost vector*/

#ifndef NAME_UBLAS_ALIAS
#define NAME_UBLAS_ALIAS
namespace ublas = boost::numeric::ublas;
#endif /*namespace alias for boost uBLAS*/

class Sphere {

Sphere(std::string, ublas::vector<float>, ublas::vector<float>, float); 
std::string getName();
ublas::vector<float> getOrigin();
ublas::vector<float> getColor(); 
float getRadius();

};

class Camera {

Camera(std::string, ublas::vector<float>, ublas::vector<float>, float, float);
std::string getName();
ublas::vector<float> getPRP();
ublas::vector<float> getVPN();
float getNearClip();
float getFarClip();

};

class Scene {

Scene(std::string, int, int, int);
std::string getName();
int getWidth();
int getHeight();
int recursionDepth();

};

class Ray {

Ray(ublas::vector<float>, ublas::vector<float>);
void setPixel(ublas::vector<float>);
void setPRP(ublas::vector<float>);
ublas::vector<float> getPixel();
ublas::vector<float> getPRP();
ublas::vector<float> rayVector(); //v = L - E, L is pixel of view plane, E is PRP
ublas::vector<float> norm();// ||v|| = sqrt( (v1)^2 + (v2)^2 +...+(vn)^2 )
ublas::vector<float> unitVector();		// U = (v/||v||)
ublas::vector<float> unitVectorScaled(float);		// sU	
ublas::vector<float> paraPos(float);// R(s) = L + sU,  L is pixel of view plane.

};

#endif  /*define OBJS_H end */
