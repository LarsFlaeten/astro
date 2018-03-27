#include "ODE.h"

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

State ODE::rates(const EphemerisTime& et, const State& s) const
{
    State sdot;
    this->operator()(s, sdot, et);
    return sdot;

}

void ODE::operator()(const State& x, State& dxdt, const EphemerisTime& et) const
{
    
    // Velocities available directly from state
    dxdt.r = x.v;

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
 


}
