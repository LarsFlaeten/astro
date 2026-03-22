#include <vector>

#include "ODE.h"
#include "Util.h"

// References:
// [1]  Spacecraft Attitude Dynamics and Control
//      Ver 3.0.0, 2008/2009, Giulio Avanzini

namespace astro {

ODE::ODE()
{}

ODE::~ODE()
{}

PosState ODE::rates(const EphemerisTime& et, const PosState& s) const
{
    PosState sdot;
    operator()(s, sdot, et);
    return sdot;
}

void ODE::operator()(const PosState& x, PosState& dxdt, const EphemerisTime& et) const
{
    dxdt.r = x.v;

    dxdt.v = Vec3(0.0);
    for (const Attractor& a : attractors)
    {
        Vec3   r = x.r - a.p;
        double R = glm::length(r);
        dxdt.v  += (-a.GM / std::pow(R, 3.0)) * r;
    }
    // TODO: add perturbations:
    // - Oblateness
    // - Atmospheric drag
    // - Solar radiation pressure
    // - Thruster forces
}

void ODE::addAttractor(const Attractor& a)
{
    attractors.push_back(a);
}

void ODE::clearAttractors()
{
    attractors.clear();
}


RotODE::RotODE(const Mat3& Ib)
    : t(0.0), t_b(0.0)
{
    setInertialMatrix(Ib);
}

RotODE::~RotODE()
{}

void RotODE::setGlobalTorque(const Vec3& _t)
{
    t = _t;
}

void RotODE::setBodyTorque(const Vec3& _tb)
{
    t_b = _tb;
}

void RotODE::setInertialMatrix(const Mat3& Ib)
{
    if (std::abs(glm::determinant(Ib)) < 1.0E-8)
        throw AstroException("Singular Matrix", "Singular matrix supplied as inertia matrix to RotODE");
    i_b     = Ib;
    i_b_inv = glm::inverse(i_b);
}


RotState RotODE::rates(const EphemerisTime& et, const RotState& rs) const
{
    RotState rs_dot;

    Quat Q     = rs.q;
    Quat Q_inv = glm::inverse(Q);
    double q0  = Q.w; // scalar part
    Vec3   q   = Vec3(Q.x, Q.y, Q.z);

    // Quaternion derivative — from [1]
    Vec3 w       = rs.w;
    double q0dot = -0.5 * glm::dot(w, q);
    Vec3   qdot  = 0.5 * (q0 * w - glm::cross(w, q));
    rs_dot.q = Quat(q0dot, qdot.x, qdot.y, qdot.z);

    // Angular velocity derivative

    // Transform global-frame torque to body frame and sum — [1] eq (56)
    Vec3 tb2  = transform(Q_inv, t, Q);
    Vec3 tbt  = t_b + tb2;

    // Angular velocity in body frame
    Vec3 wb = transform(Q_inv, w, Q);

    // [1] eq (4): L_dot = tau - w × (I * w)
    // TODO: For time-varying I, subtract I_dot * wb
    Vec3 Lb_dot = tbt - glm::cross(wb, i_b * wb);
    Vec3 wbdot  = i_b_inv * Lb_dot;
    rs_dot.w    = transform(Q, wbdot, Q_inv);

    return rs_dot;
}

} // namespace astro
