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
    // Pi function, and a static variable (no need to evaluate pi every time..)
    double  pi();
    const double PI = pi();



	// Degrees / radians conversion
	double  degreesPerRadian();
    double  radiansPerDegree();
    const double DEGPERRAD = degreesPerRadian();
    const double RADPERDEG = radiansPerDegree();    

    // convert to a position from Right Ascension, Declination
    void    raDecToVec(double range, double ra, double dec, ork::vec3d* res);

    // convert from a position to Right Ascension, Declination
    void    vecToRaDec(const ork::vec3d& pos, double* range, double* ra, double* dec);

}



#endif
