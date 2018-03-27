#ifndef _ASTRO_RKF45_H_
#define _ASTRO_RKF45_H_

#include <vector>
#include "Time.h"
#include "State.h"
#include "ODE.h"
namespace astro
{

class RKF45
{
public:
    struct Result
    {
        State s;
        EphemerisTime et;
        TimeDelta dt_next;
    };


    static Result doStep(const ODE& ode, const State& s, const EphemerisTime& et, const TimeDelta& dt);            

    static std::vector<Result> doSteps(const ODE& ode, const State& s, const EphemerisTime& et0, const EphemerisTime& et1, const TimeDelta& dt);            



    static void setTolerance(double _tol);

private:
    // The tolerance to be used when estimating next time step
    // Default is 1.0E-8;
    static double tol;

    static const std::vector<double> a;
    static const std::vector< std::vector<double> > b;
    static const std::vector<double> c4;
    static const std::vector<double> c5; 

    static const double eps;
	static const double h_min;
};




}


#endif
