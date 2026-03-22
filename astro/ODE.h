#ifndef _ASTRO_ODE_H_
#define _ASTRO_ODE_H_

#include "State.h"
#include "Time.h"
#include "Exceptions.h"
#include <vector>

namespace astro {


class Attractor
{
public:
    Vec3   p;  // Position relative to the frame of reference
    double GM;
};

// Differential equation for translation. Separated from rotation because:
// - They are largely uncoupled
// - They have very different optimal time steps
// - This allows different propagation regimes (e.g. numerical for rotation,
//   Keplerian for translation)
class ODE
{
public:
    ODE();
    virtual ~ODE();

    // Evaluates force terms and returns derivatives of the state vector
    PosState rates(const EphemerisTime& et, const PosState& s) const;

    // Callable interface required by ODE solvers
    virtual void operator()(const PosState& x, PosState& dxdt, const EphemerisTime& et) const;

    virtual void addAttractor(const Attractor& a);
    virtual void clearAttractors();

private:
    std::vector<Attractor> attractors;
};


class RotODE
{
public:
    // Ib is the inertia matrix in body frame. Default is the identity matrix.
    explicit RotODE(const Mat3& Ib = Mat3(1.0));
    virtual ~RotODE();

    // Force function for rotational dynamics
    RotState rates(const EphemerisTime& et, const RotState& s) const;

    // Set global-frame torque (e.g. gravity gradient)
    void setGlobalTorque(const Vec3& t);

    // Set body-frame torque (e.g. thrusters)
    void setBodyTorque(const Vec3& tb);

    // Set the inertia matrix in body frame
    void setInertialMatrix(const Mat3& Ib);

private:
    Vec3 t;       // global torque
    Vec3 t_b;     // body-frame torque
    Mat3 i_b;     // body-frame inertia matrix
    Mat3 i_b_inv; // inverse inertia matrix
};


} // namespace astro

#endif
