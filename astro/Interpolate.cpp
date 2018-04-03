#include "Interpolate.h"

namespace astro
{

void   hermite(const PosState& s1, const EphemerisTime& et1, const PosState& s2, const EphemerisTime& et2, const EphemerisTime& etx, PosState& out)
{
    // interpolation on arbitrary interval:
    // https://en.wikipedia.org/wiki/Cubic_Hermite_spline#Interpolation_on_an_arbitrary_interval
    double span = (et2 - et1).value;
    double t = (etx - et1).value / span;

    // Hermite basis functions
    double h00 =  2.0*pow(t,3.0) - 3.0*pow(t,2.0) + 1;
    double h01 = -2.0*pow(t,3.0) + 3.0*pow(t,2.0);
    double h10 =      pow(t,3.0) - 2.0*pow(t,2.0) + t;
    double h11 =      pow(t,3.0) -     pow(t,2.0);

   
    out.r = h00*s1.r + h10*span*s1.v + h01*s2.r + h11*span*s2.v;
    
    // Derivatives:
    // Ref discussion here:
    // https://math.stackexchange.com/questions/2444650/cubic-hermite-spline-derivative
    h00 =  6.0*pow(t,2.0) - 6.0*t;
    h01 = -6.0*pow(t,2.0) + 6.0*t;
    h10 =  3.0*pow(t,2.0) - 4.0*t + 1;
    h11 =  3.0*pow(t,2.0) - 2.0*t;
    out.v = (h00*s1.r + h10*span*s1.v + h01*s2.r + h11*span*s2.v)/span;

}    





}
