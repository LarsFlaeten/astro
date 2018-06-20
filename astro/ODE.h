#ifndef _ASTRO_ODE_H_
#define _ASTRO_ODE_H_

#include "State.h"
#include "Time.h"
#include "Exceptions.h"

#include <ork/math/mat3.h>

namespace astro
{


class Attractor
{
public:
    ork::vec3d  p; // Position relative to the frame of reference used
    double      GM;
};
// This is the differential equation for translation. Translation is separated
// from rotation as these will be integrated separately:
// - They are not coupled or at least lightly coupled
// - They wil have very different optimal time steps. E.g. an hyperbolic escape
//   trajectory from earth with RKF78 will have dt ~ 1.0E6 / 1.0E7 @ SOI, while
//   rotation may need to be integrated with dt many orders below this
// - This also allows different propagation regimes, for example numerical integration
//   for rotations, while using Kepler/conincs for translation.
class ODE
{
public:
    ODE();
    virtual ~ODE();

 

    // The force function which evaluates the given force terms
    // and returns the derivatives of the state vector given
    PosState rates(const EphemerisTime& et, const PosState& s) const;

    // Essentially same method as rates, but odeint needs it slightly different
    virtual void operator()(const PosState& x, PosState& dxdt, const EphemerisTime& et) const;

    virtual void addAttractor(const Attractor& a);
    virtual void clearAttractors();

private:
    std::vector<Attractor> attractors;
};

class RotODE
{
public:
    // CTOR
    // Ib is the inertial matrix
    // Default is the identity matrix
    RotODE(const ork::mat3d& Ib = ork::mat3d::IDENTITY);

    virtual ~RotODE();

    // The force function for rotations
    RotState rates(const EphemerisTime& et, const RotState& s) const;
    
    // Assume constant torque over the time step
    // TODO: Implement linear varying torque over the time step?
    // Set global frame torque
    // (Gravity gradient etc)
    void setGlobalTorque(const ork::vec3d& t);

    // Set body frame torque
    // (From thrusters etc)
    void setBodyTorque(const ork::vec3d& tb);

    // Set the inertial matrix in body frame
    void setInertialMatrix(const ork::mat3d& Ib);

private:

    vec3d   t;      // global torque
    vec3d   t_b;    // body frame torque
    mat3d   i_b;    // body fram inertal matrix
    mat3d   i_b_inv;// inverse inertial matrix

};


}


#endif
