#ifndef _ASTRO_INTERPOLATE_H_
#define _ASTRO_INTERPOLATE_H_

#include "State.h"

namespace astro
{


// Performs a Hermite spline interpolation between state 1 and 2,
// giving continous first order derivatives (velocities)
// f is a fraction in "time", betwen 0 (s1) and 1 (s2)
void   hermite(const PosState& s1, const PosState& s2, double f, PosState& out);


}


#endif

