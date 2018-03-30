#include "RKF78.h"
#include "Exceptions.h"

#include <limits>
#include <algorithm>
#include <iostream>
#include <sstream>

namespace astro
{

double RKF78::tol = 1.0E-8;

int RKF78::allowable_tries = 10;

rkf78_controlled_stepper RKF78::stepper = boost::numeric::odeint::make_controlled(RKF78::tol, RKF78::tol, rkf78_error_stepper() );



void    RKF78::setTolerance(double _tol)
{
    if(_tol < 0.0)
        throw AstroException("Zero or negative tolerance not allowed for RK methods");
    RKF78::tol = _tol;

    // We have to make a new stepper with updated tolerance, as there is no way 
    // in odeint to update the tolerances of a stepper instance..
    RKF78::stepper = boost::numeric::odeint::make_controlled(RKF78::tol, RKF78::tol, rkf78_error_stepper ());

}

RKF78::Result RKF78::doStep(const ODE& ode, const State& s, const EphemerisTime& et, const TimeDelta& dt)
{
    double et_d = et.getETValue(); 
    double DT = dt.value;
    int tries = 0;
    State si = s;

    while(stepper.try_step(ode, si, et_d, DT) == boost::numeric::odeint::controlled_step_result::fail)
    {
        tries++;
        if(tries>allowable_tries)
            throw AstroException("RKF78 did not converge within allowable tries");

    }
    EphemerisTime et1(et_d);
    TimeDelta dt_next(DT);

    // prepare result:
    Result res;
    res.dt_next = dt_next;
    res.et = et1;
    res.s = si;
    res.numTries = tries;
    
	return res;

}
 
std::vector<RKF78::Result> RKF78::doSteps(const ODE& ode, const State& s, const EphemerisTime& et0, const EphemerisTime& et1, const TimeDelta& dt)
{
    std::vector<Result> res;
    // Initial value:
    res.push_back( RKF78::Result({s, et0, dt}));


    while(res.back().et < et1)
    {
        res.push_back(doStep(ode, res.back().s, res.back().et, res.back().dt_next));
        
        if(res.back().et + res.back().dt_next > et1)
            res.back().dt_next = et1 - res.back().et;
         
    }


    return std::move(res);
}            


}
