#include "ODE.h"

// References:
// [1]  Spacecraft Attitude Dynamics and Control
//      Ver 3.0.0, 2008/2009 Giulio Avanzini

namespace astro
{

ODE::ODE(double _mu)
    : mu(_mu)
{

}

ODE::~ODE()
{

}

void ODE::setMu(double _mu)
{
    mu = _mu;
}   

PosState ODE::rates(const EphemerisTime& et, const PosState& s) const
{
    PosState sdot;
    this->operator()(s, sdot, et);
    return sdot;

}

void ODE::operator()(const PosState& x, PosState& dxdt, const EphemerisTime& et) const
{
    // Set r_dot:
    // Velocities available directly from state
    dxdt.r = x.v;

    // Set v_dot:
    // Main gravitational force
    double r = x.r.length();    
    dxdt.v = x.r/pow(r, 3.0);
    dxdt.v *= -mu;

    // TODO:add other things here:
    // -Oblateness effect
    // -Atmospheric drag
    // -solar radiation drag
    // -thruster forces


}
 
RotODE::RotODE()
{

}

RotODE::~RotODE()
{

}

RotState RotODE::rates(const EphemerisTime& et, const RotState& rs) const
{
    RotState rs_dot;

    quatd Q = rs.q;
    double q0 = Q.w; // w is the scalar in ork/quatd
    vec3d q = vec3d(Q.x, Q.y, Q.z);

    // Derivative of the quaternion
    // From [1]
    vec3d w = rs.w;
    double q0dot = -0.5*w.dotproduct(q);
    vec3d qdot = 0.5*(q0*w - w.crossProduct(q));
    quatd Qdot = quatd(qdot.x, qdot.y, qdot.z, q0dot);
    rs_dot.q = Qdot;

    // Derivative of w
    // From torques etc
    // TODO: Implement. We then need inertial tensor and setters
    // for torque (body fixed and global)
    //rs.dot.w = I_inv*(torque +++ ) or somehting
    rs_dot.w = vec3d::ZERO;
    

    return rs_dot;
}

}
