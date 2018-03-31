#ifndef _ASTRO_PROPAGATOR_H_
#define _ASTRO_PROPAGATOR_H_

#include "State.h"
#include "Time.h"
#include "ODE.h"
#include "RKF45.h"
#include "RKF78.h"
#include "RK1_4.h"
namespace astro
{

struct SimpleResult
{
    PosState s;
    TimeDelta dt_next;
};



template< typename ODEType, typename Solver, typename Result = typename Solver::Result >
class Propagator
{
public:
    Propagator(const ODEType& ode);
    ~Propagator();
   
    // Perform one numerical integrarion step
    // State s - State to be integrated
    // et - time at which the state is given
    // dt - stepsize, is often given by a previous integration result 
    Result doStep(const PosState& s, const EphemerisTime& et, const TimeDelta& dt);

    // Perform a step sequence from et0 to et1:
    // State s - State to be integrated
    // et0 - time at which the initial state is given
    // et1 - time at the end of integration
    // dt - initial stepsize 
    std::vector<Result> doSteps(const PosState& s, const EphemerisTime& et0, const EphemerisTime& et1, const TimeDelta& dt);
    

    Solver& getSolver();

private:
    const ODEType& ode;

};

template< typename ODEType, typename Solver, typename Result >
Propagator<ODEType, Solver, Result>::Propagator(const ODEType& ode)
    : ode(ode)
{

}

template< typename ODEType, typename Solver, typename Result >
Propagator<ODEType, Solver, Result>::~Propagator()
{

}



template< typename ODEType, typename Solver, typename Result >
Result Propagator<ODEType, Solver, Result>::doStep(const PosState& s, const EphemerisTime& et, const TimeDelta& dt)
{
    return Solver::doStep(ode, s, et, dt);
}
 
template< typename ODEType, typename Solver, typename Result >
std::vector<Result> Propagator<ODEType, Solver, Result>::doSteps(const PosState& s, const EphemerisTime& et0, const EphemerisTime& et1, const TimeDelta& dt)
{
    return std::move(Solver::doSteps(ode, s, et0, et1, dt));
}
 
 
}



#endif

