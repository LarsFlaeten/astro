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
    RotState rsn_dot = rode.rates(et, rs);

    //quatd qn_dot = q_wn*qn;
    //qn_dot.x *= 0.5;
    //qn_dot.y *= 0.5;
    //qn_dot.z *= 0.5;
    //qn_dot.w *= 0.5;

    //std::cout << qn_dot << std::endl;

    // (57)
    // angular velocities at n+1/4 and n+1/2
    vec3d wndot = rsn_dot.w;

    // Transform wndot to body coordinates wbndot
    vec3d wbndot = transform(qn_inv, wndot, qn);

    vec3d wbn14 = wbn + 0.25*wbndot*DT;
    vec3d wbn12 = wbn + 0.5 *wbndot*DT;

    // (58)
    // Angular velocity in global frame at n+1/4
    vec3d wn14 = transform(qn, wbn14, qn_inv);

    // (59)
    // Calculate predicted q'n12:
    double wn14_l = wn14.length();
    double F = wn14_l*DT*0.25;
    vec3d tmp = wn14*(1.0/wn14_l);
    tmp *= sin(F);
    // Add guard against zero angular velocity:
    if(wn14_l < 1.0E-8)
    {
        F = 0.0;
        tmp = vec3d::ZERO;
    }
    quatd qprime_n12 = quatd(tmp.x, tmp.y, tmp.z, cos(F))*qn;
    
    // (60)
    vec3d wn12 = transform(qprime_n12, wbn12, qprime_n12.inverse());

    // Intermdiate step, we can now get derivates for n+1/2:
    // Note: With current implementation of RotODE, the torque is
    // constant and hence w_dot is roughly constant over the time step
    // This may change in the future by applying other RotODEs.
    RotState rprime_n12(qprime_n12, wn12);
    TimeDelta dt12(0.5*DT);
    RotState rsn12_dot = rode.rates(et+dt12, rprime_n12);

    // (61)
    // Calculate q_n+1
    double wn12_l = wn12.length();
    F = wn12_l*DT*0.5;
    tmp = wn12*(1.0/wn12_l);
    tmp *= sin(F);
    // Add guard against zero velocity
    if(wn12_l < 1.0E-8)
    {
        F = 0.0;
        tmp = vec3d::ZERO;
    }
    quatd qn1 = quatd(tmp.x, tmp.y, tmp.z, cos(F))*qn;
    quatd qn1_inv = qn1.inverse();

    // (62)
    // Angular velocities at n+1:
    vec3d wbn12dot = rsn12_dot.w;
    vec3d wbn1 = wbn + wbn12dot*DT;
    // (63)
    // Transform to global frame:
    vec3d wn1 = transform(qn1, wbn1, qn1_inv);

    Result res;
    res.rs.q = qn1;
	res.rs.w = wn1;
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
