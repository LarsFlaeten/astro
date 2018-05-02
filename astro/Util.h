#ifndef _ASTRO_UTIL_H
#define _ASTRO_UTIL_H

// Utility functions used by many other parts of the library

// We use ork::math
// Need to have ORK_API defined. This hsould be covered by the overlying application, such as ssim?
#ifndef ORK_API
#define ORK_API
#endif
#include <ork/math/vec3.h>
#include <ork/math/quat.h>



#include <sstream>


namespace astro
{
    // Pi function, and a few static variables (no need to evaluate pi every time..)
    double  pi();
    const double PI = pi();
    const double TWOPI = 2.0*pi();
    const double PIHALF = pi()*0.5;


	// Degrees / radians conversion
	double  degreesPerRadian();
    double  radiansPerDegree();
    const double DEGPERRAD = degreesPerRadian();
    const double RADPERDEG = radiansPerDegree();    

    // convert to a position from Right Ascension, Declination
    void    raDecToVec(double range, double ra, double dec, ork::vec3d* res);

    // convert from a position to Right Ascension, Declination
    void    vecToRaDec(const ork::vec3d& pos, double* range, double* ra, double* dec);

	//Normalizes any number to an arbitrary range 
    double wrap( const double value, const double start, const double end ); 
    
    // Transforms a vector by a quaternion (and its inverse)
    // Usage (if Q denotes the orientation of a body):
    // v_b = transform(Q_inv, v,   Q)
    // v   = transform(Q,     v_b, Q_inv)
    ork::vec3d transform(const ork::quatd& q1, const ork::vec3d& v, const ork::quatd& q2);


    // Prettyprinters 
    std::ostream& operator << (std::ostream& os, const ork::vec3d& v);
    std::ostream& operator << (std::ostream& os, const ork::quatd& q);

}





#endif
