#include "PCDM.h"
#include "Util.h"

namespace astro
{

PCDM::Result PCDM::doStep(const RotODE& rode, const RotState& rs, const EphemerisTime& et, const TimeDelta& dt)
{
    double et_d = et.getETValue(); 
    double DT = dt.value;
    
    // Rotation at time n
    quatd qn = rs.q;
    quatd qn_inv = qn.inverse();

    // Angular velocity at time n, quaternized
    quatd q_wn = quatd(rs.w.x, rs.w.y, rs.w.z, 0);

    // Transform angular velocities (and later any global torques) to body fixed:
    quatd q_wbn = qn_inv*q_wn*qn;
    vec3d wbn = vec3d(q_wbn.x, q_wbn.y, q_wbn.z); 
    
    // Get the derivatives for time/step n:
    // TODO: Drop the ODe and calculate directly here?
    RotState rsn_dot = rode.rates(et, rs);
    //std::cout << rsn_dot.q << std::endl;

    quatd qn_dot = q_wn*qn;
    qn_dot.x *= 0.5;
    qn_dot.y *= 0.5;
    qn_dot.z *= 0.5;
    qn_dot.w *= 0.5;

    //std::cout << qn_dot << std::endl;

    // (57)
    // angular velocities at n+1/4 and n+1/2
    vec3d wndot = rsn_dot.w;
    vec3d wbn14 = wbn + 0.25*wndot*DT;
    vec3d wbn12 = wbn + 0.5 *wndot*DT;

    // (58)
    // Angular velocity in global frame at n+1/4
    quatd q_wbn14 = quatd(wbn14.x, wbn14.y, wbn14.z, 0);
    quatd q_wn14 = qn*q_wbn14*qn_inv;
    vec3d wn14 = vec3d(q_wn14.x, q_wn14.y, q_wn14.z);

    // (59)
    // Calculate predicted q'n12:
    double wn14_l = wn14.length();
    double F = wn14_l*DT*0.25;
    vec3d tmp = wn14*(1.0/wn14_l);
    tmp *= sin(F);
    quatd qprime_n12 = quatd(tmp.x, tmp.y, tmp.z, cos(F))*qn;
    
    // (60)
    quatd q_wbn12 = quatd(wbn12.x, wbn12.y, wbn12.z, 0);
    quatd q_wn12 = qprime_n12 * q_wbn12 * qprime_n12.inverse();
    vec3d wn12 = vec3d(q_wn12.x, q_wn12.y, q_wn12.z);

    // (61)
    // Calculate q_n+1
    double wn12_l = wn12.length();
    F = wn12_l*DT*0.5;
    tmp = wn12*(1.0/wn12_l);
    tmp *= sin(F);
    quatd qn1 = quatd(tmp.x, tmp.y, tmp.z, cos(F))*qn;


    // (62)
    // Angular velocities at n+1:
    // TODO: we need wbn12_dot!!
    //vec3d wbn1 = wbn + w

    //std::cout << rsn_dot.q.x << ", " << rsn_dot.q.y << ", " << rsn_dot.q.z << ", " << rsn_dot.q.w << std::endl;
    //std::cout << qn_dot.x << ", " << qn_dot.y << ", " << qn_dot.z << ", " << qn_dot.w << std::endl;
    
    



    Result res;
    res.rs.q = qn1;
	res.rs.w = rs.w; // TODO: get from above
    res.et = et + dt;

    return res;

}
 
std::vector<PCDM::Result> PCDM::doSteps(const RotODE& rode, const RotState& rs, const EphemerisTime& et0, const EphemerisTime& et1, const TimeDelta& dt)
{
    std::vector<Result> res;
    // Initial value:
    res.push_back( PCDM::Result({rs, et0}));


    TimeDelta dti = dt;

    while(res.back().et < et1)
    {
        res.push_back(doStep(rode, res.back().rs, res.back().et, dti));
        
        if(res.back().et + dti > et1)
            dti = et1 - res.back().et;
         
    }


    return std::move(res);
}            


}
