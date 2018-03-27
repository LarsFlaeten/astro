#include "Interpolate.h"

namespace astro
{

void   hermite(const State& s1, const State& s2, double f, State& out)
{
    double h1 =  2.0*pow(f,3.0) - 3.0*pow(f,2.0) + 1;
    double h2 = -2.0*pow(f,3.0) + 3.0*pow(f,2.0);
    double h3 =      pow(f,3.0) - 2.0*pow(f,2.0) + f;
    double h4 =      pow(f,3.0) -     pow(f,2.0);

   
    out.r = s1.r*h1 + s2.r*h2 + s1.v*h3 + s2.v*h4;
    
    // Derivatives:
    h1 =  6.0*pow(f,2.0) - 6.0*f;
    h2 = -6.0*pow(f,2.0) + 6.0*f;
    h3 =  3.0*pow(f,2.0) - 4.0*f + 1;
    h4 =  3.0*pow(f,2.0) - 2.0*f;
    out.v = s1.r*h1 + s2.r*h2 + s1.v*h3 + s2.v*h4;

    



}

}
