#ifndef _ASTRO_RKF78_H_
#define _ASTRO_RKF78_H_

#include <vector>
#include "Time.h"
#include "State.h"
#include "ODE.h"

namespace astro {

// Runge-Kutta-Fehlberg 7(8) adaptive integrator.
// Butcher tableau from Fehlberg (1969), 13 stages.
// Uses the 8th-order solution for propagation and the difference
// between 7th and 8th order as the error estimate for step-size control.
class RKF78
{
public:
    struct Result
    {
        PosState     s;
        EphemerisTime et;
        TimeDelta    dt_next;
        int          numTries;
    };

    static Result doStep(const ODE& ode, const PosState& s, const EphemerisTime& et, const TimeDelta& dt);

    static std::vector<Result> doSteps(const ODE& ode, const PosState& s,
                                       const EphemerisTime& et0, const EphemerisTime& et1,
                                       const TimeDelta& dt);

    static void setTolerance(double tol);

private:
    static double tol;

    // Butcher tableau — nodes, coupling matrix, 7th and 8th order weights
    static const double              c[13];
    static const std::vector<double> b[13]; // b[i] has i entries (lower triangle)
    static const double              ch7[13]; // 7th order weights
    static const double              ch8[13]; // 8th order weights

    static const double eps;
};

} // namespace astro

#endif
