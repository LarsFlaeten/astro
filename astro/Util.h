#ifndef _ASTRO_UTIL_H
#define _ASTRO_UTIL_H

// Utility functions used by many other parts of the library

// We use ork::math
// Need to have ORK_API defined. This hsould be covered by the overlying application, such as ssim?
#ifndef ORK_API
#define ORK_API
#endif
#include <ork/math/vec3.h>

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
 
}



#endif
