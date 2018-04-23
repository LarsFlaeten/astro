#ifndef _ASTRO_ODE_H_
#define _ASTRO_ODE_H_

#include "State.h"
#include "Time.h"

namespace astro
{

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
    ODE(double _mu);
    virtual ~ODE();

    // Set the current gravitational parameter to be used
    void    setMu(double _mu);
 

    // The force function which evaluates the given force terms
    // and returns the derivatives of the state vector given
    PosState rates(const EphemerisTime& et, const PosState& s) const;

    // Essentially same method as rates, but odeint needs it slightly different
    virtual void operator()(const PosState& x, PosState& dxdt, const EphemerisTime& et) const;

private:
    double mu;
};

class RotODE
{
public:
    RotODE();
    virtual ~RotODE();

    // The force function for rotations
    RotState rates(const EphemerisTime& et, const RotState& s) const;

};


}


#endif
