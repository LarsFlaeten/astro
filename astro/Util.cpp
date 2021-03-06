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

    void    raDecToVec(double range, double ra, double dec, mork::vec3d* res)
    {
        radrec_c(range, ra, dec, &(res->x));
    }
      
    void    vecToRaDec(const mork::vec3d& pos, double* range, double* ra, double* dec)
    {
        recrad_c(&(pos.x), range, ra, dec);
    }

    //Normalizes any number to an arbitrary range 
    //by assuming the range wraps around when going below min or above max 
    double wrap( const double value, const double start, const double end ) 
    {
        const double width       = end - start   ;   // 
        const double offsetValue = value - start ;   // value relative to 0

        return ( offsetValue - ( floor( offsetValue / width ) * width ) ) + start ;
        // + start to reset back to start of original range
    }

    mork::vec3d transform(const mork::quatd& q1, const mork::vec3d& v, const mork::quatd& q2)
    {
        mork::quatd q_v = mork::quatd(v.x, v.y, v.z, 0.0);
        mork::quatd res_q_v = q1 * q_v * q2;
        return mork::vec3d(res_q_v.x, res_q_v.y, res_q_v.z);
    }



    std::ostream& operator << (std::ostream& os, const mork::vec3d& v)
    {
        os << "[ " << v.x << ", " << v.y << ", " << v.z << " ]";
        return os;
    }

    std::ostream& operator << (std::ostream& os, const mork::quatd& q)
    {
        os << "[ " << q.x << "i + " << q.y << "j + " << q.z << "k + " << q.w << " ]";
        return os;
    }
    
    std::ostream& operator << (std::ostream& os, const mork::mat3d& m)
    {
        for(int i = 0; i < 3; ++i)
            os << "[ " << m[i][0] << ", " << m[i][1] << ", " << m[i][2] << " ]\n";
        return os;
    } 



 
}

