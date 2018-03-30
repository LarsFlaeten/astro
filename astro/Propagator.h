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



template< typename Solver, typename Result >
class Propagator
{
public:
    Propagator(const ODE& ode);
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
    const ODE& ode;

};

template< typename Solver, typename Result >
Propagator<Solver,Result>::Propagator(const ODE& ode)
    : ode(ode)
{

}

template< typename Solver, typename Result >
Propagator<Solver,Result>::~Propagator()
{

}



template< typename Solver, typename Result >
Result Propagator<Solver,Result>::doStep(const PosState& s, const EphemerisTime& et, const TimeDelta& dt)
{
    return Solver::doStep(ode, s, et, dt);
}
 
template< typename Solver, typename Result >
std::vector<Result> Propagator<Solver,Result>::doSteps(const PosState& s, const EphemerisTime& et0, const EphemerisTime& et1, const TimeDelta& dt)
{
    return std::move(Solver::doSteps(ode, s, et0, et1, dt));
}
 
 
}



#endif

