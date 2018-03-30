#ifndef _ASTRO_RKF78_H_
#define _ASTRO_RKF78_H_

#include <vector>
#include "Time.h"
#include "State.h"
#include "ODE.h"

#include <boost/numeric/odeint.hpp>

using namespace boost::numeric::odeint;

namespace astro
{

typedef boost::numeric::odeint::runge_kutta_fehlberg78< astro::PosState, double, astro::PosState, double, boost::numeric::odeint::vector_space_algebra > rkf78_error_stepper;

typedef boost::numeric::odeint::controlled_runge_kutta < rkf78_error_stepper > rkf78_controlled_stepper;

class RKF78
{
public:
    struct Result
    {
        PosState s;
        EphemerisTime et;
        TimeDelta dt_next;
        int numTries;
    };


    static Result doStep(const ODE& ode, const PosState& s, const EphemerisTime& et, const TimeDelta& dt);            

    static std::vector<Result> doSteps(const ODE& ode, const PosState& s, const EphemerisTime& et0, const EphemerisTime& et1, const TimeDelta& dt);            



    static void setTolerance(double _tol);

private:
    // The tolerance to be used when estimating next time step
    // Default is 1.0E-8;
    static double tol;
    static int allowable_tries;
    static rkf78_controlled_stepper stepper;
};




}


#endif
