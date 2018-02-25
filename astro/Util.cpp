#include "Util.h"

// Wrap a lot of the Spice utility functions
#include <cspice/SpiceUsr.h>

namespace astro
{
    double  pi()
    {
        return pi_c();
    }

    double  degreesPerRadian()
    {
        return   dpr_c();
    }
    
    double  radiansPerDegree()
    {
        return 1.0 / dpr_c();
    }

    void    raDecToVec(double range, double ra, double dec, ork::vec3d* res)
    {
        radrec_c(range, ra, dec, &(res->x));
    }
      
    void    vecToRaDec(const ork::vec3d& pos, double* range, double* ra, double* dec)
    {
        recrad_c(&(pos.x), range, ra, dec);
    }

}
