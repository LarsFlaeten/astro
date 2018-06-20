#include "ODE.h"
#include "Util.h"

// References:
// [1]  Spacecraft Attitude Dynamics and Control
//      Ver 3.0.0, 2008/2009 Giulio Avanzini

namespace astro
{

ODE::ODE()
{

}

ODE::~ODE()
{

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
    // Main gravitational force from attractors
    dxdt.v = ork::vec3d::ZERO;
    for(const Attractor& a : attractors) {
        vec3d ac = ork::vec3d::ZERO;
        vec3d r = x.r - a.p;
        double R = r.length();    
        ac = r/pow(R, 3.0);
        ac *= -a.GM;
        dxdt.v += ac;
    }
    // TODO:add other things here:
    // -Oblateness effect
    // -Atmospheric drag
    // -solar radiation drag
    // -thruster forces


}

void ODE::addAttractor(const Attractor& a)
{
    attractors.push_back(a);
}

void ODE::clearAttractors()
{
    attractors.clear();
}


RotODE::RotODE(const mat3d& Ib)
{
    setInertialMatrix(Ib);
    t = ork::vec3d::ZERO;
    t_b = ork::vec3d::ZERO; 

}

RotODE::~RotODE()
{

}

void RotODE::setGlobalTorque(const ork::vec3d& _t)
{
    t = _t;
}


void RotODE::setBodyTorque(const ork::vec3d& _tb)
{
    t_b = _tb;
}

void RotODE::setInertialMatrix(const ork::mat3d& Ib)
{
    if(fabs(Ib.determinant()) < 1.0E-8)
        throw AstroException("Singular Matrix", "Singular Matrix supplied as inertial matirx to RotODE");

    i_b = Ib;
    i_b_inv = i_b.inverse();

}


RotState RotODE::rates(const EphemerisTime& et, const RotState& rs) const
{
    RotState rs_dot;

    quatd Q = rs.q;
    quatd Q_inv = Q.inverse();
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
    
    // Transform global frame torque to body frame, and add the two:
    vec3d tb2 = transform(Q_inv, t, Q); // [1] (56)
    vec3d tbt = t_b + tb2; // Total body torque

    // Get w in body frame, wb:
    vec3d wb = transform(Q_inv, w, Q);
    

    // [1] (4)
    vec3d Lb_dot = tbt - wb.crossProduct(i_b*wb); 
    // TODO: For varying I we also need to substract Ib_dot*wb
    vec3d wbdot = i_b_inv * (Lb_dot /* - i_b_dot*wb*/);
    vec3d wdot = transform(Q, wbdot, Q_inv);

    rs_dot.w = wdot;
    

    return rs_dot;
}

}
