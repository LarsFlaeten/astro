#ifndef _ASTRO_ODE_H_
#define _ASTRO_ODE_H_

#include "State.h"
#include "Time.h"

namespace astro
{


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

}


#endif
