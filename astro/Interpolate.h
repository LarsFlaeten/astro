#ifndef _ASTRO_INTERPOLATE_H_
#define _ASTRO_INTERPOLATE_H_

#include "State.h"
#include "Time.h"

namespace astro
{


// Performs a Hermite spline interpolation between state 1 and 2,
// giving continous first order derivatives (velocities)
// s1  - State at t1;
// et1 - Time at t1;
// s2  - State at t2;
// et2 - Time at t2
// etx - Time at which to fint interpolated state; et1 < etx < et2
// out - The interpolated state
void   hermite(const PosState& s1, const EphemerisTime& et1, const PosState& s2, const EphemerisTime& et2, const EphemerisTime& etx, PosState& out);


}


#endif

