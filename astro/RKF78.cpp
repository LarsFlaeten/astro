#include "RKF78.h"
#include "Exceptions.h"

#include <limits>
#include <cmath>
#include <algorithm>
#include <sstream>

// Fehlberg 7(8) Butcher tableau — 13 stages.
// Reference: Fehlberg, E. (1969), NASA TR R-287
// Coefficients verified against boost::numeric::odeint::runge_kutta_fehlberg78.

namespace astro {

// ── static data ─────────────────────────────────────────────────────────────

double RKF78::tol = 1.0E-8;

const double RKF78::eps = std::numeric_limits<double>::epsilon();

// Nodes
const double RKF78::c[13] = {
    0,
    2.0/27.0,
    1.0/9.0,
    1.0/6.0,
    5.0/12.0,
    1.0/2.0,
    5.0/6.0,
    1.0/6.0,
    2.0/3.0,
    1.0/3.0,
    1.0,
    0.0,
    1.0
};

// Coupling coefficients (lower-triangular Butcher matrix).
// b[i] holds the i non-zero entries for stage i+1.
const std::vector<double> RKF78::b[13] = {
    {},                                                                                      // k1
    { 2.0/27.0 },                                                                            // k2
    { 1.0/36.0,   1.0/12.0 },                                                               // k3
    { 1.0/24.0,   0.0,            1.0/8.0 },                                                // k4
    { 5.0/12.0,   0.0,           -25.0/16.0,    25.0/16.0 },                               // k5
    { 1.0/20.0,   0.0,            0.0,           1.0/4.0,       1.0/5.0 },                 // k6
    {-25.0/108.0, 0.0,            0.0,          125.0/108.0,  -65.0/27.0,   125.0/54.0 }, // k7
    { 31.0/300.0, 0.0,            0.0,           0.0,           61.0/225.0,  -2.0/9.0,
       13.0/900.0 },                                                                         // k8
    {  2.0,       0.0,            0.0,          -53.0/6.0,    704.0/45.0,  -107.0/9.0,
       67.0/90.0,  3.0 },                                                                   // k9
    {-91.0/108.0, 0.0,            0.0,           23.0/108.0, -976.0/135.0,  311.0/54.0,
      -19.0/60.0,  17.0/6.0,     -1.0/12.0 },                                              // k10
    {2383.0/4100.0, 0.0,          0.0,          -341.0/164.0, 4496.0/1025.0,-301.0/82.0,
     2133.0/4100.0, 45.0/82.0,   45.0/164.0,   18.0/41.0 },                               // k11
    {   3.0/205.0, 0.0,           0.0,           0.0,          0.0,          -6.0/41.0,
       -3.0/205.0,-3.0/41.0,      3.0/41.0,     6.0/41.0,     0.0 },                      // k12
    {-1777.0/4100.0, 0.0,         0.0,          -341.0/164.0, 4496.0/1025.0,-289.0/82.0,
     2193.0/4100.0,  51.0/82.0,   33.0/164.0,   12.0/41.0,   0.0,           1.0 }         // k13
};

// 7th-order weights
const double RKF78::ch7[13] = {
    41.0/840.0, 0, 0, 0, 0,
    34.0/105.0,
     9.0/35.0,
     9.0/35.0,
     9.0/280.0,
     9.0/280.0,
    41.0/840.0,
    0, 0
};

// 8th-order weights
const double RKF78::ch8[13] = {
    0, 0, 0, 0, 0,
    34.0/105.0,
     9.0/35.0,
     9.0/35.0,
     9.0/280.0,
     9.0/280.0,
    0,
    41.0/840.0,
    41.0/840.0
};

// ── helpers ──────────────────────────────────────────────────────────────────

static PosState stageState(const PosState& s, double h,
                           const std::vector<double>& brow,
                           const std::vector<PosState>& k)
{
    PosState result = s;
    for (size_t j = 0; j < brow.size(); ++j)
        result += k[j] * (h * brow[j]);
    return result;
}

// ── interface ────────────────────────────────────────────────────────────────

void RKF78::setTolerance(double _tol)
{
    if (_tol <= 0.0)
        throw AstroException("Zero or negative tolerance not allowed for RK methods");
    tol = _tol;
}

RKF78::Result RKF78::doStep(const ODE& ode, const PosState& s,
                             const EphemerisTime& et, const TimeDelta& dt)
{
    const double h  = dt.value;
    const double ti = et.getETValue();

    // Evaluate 13 stage derivatives
    std::vector<PosState> k(13);
    for (int i = 0; i < 13; ++i)
    {
        PosState si = stageState(s, h, b[i], k);
        k[i] = ode.rates(EphemerisTime(ti + c[i] * h), si);
    }

    // Error estimate: h * (ch7 - ch8) = h * 41/840 * (k[0]+k[10] - k[11]-k[12])
    PosState te;
    te.r = Vec3(0.0);
    te.v = Vec3(0.0);
    for (int i = 0; i < 13; ++i)
    {
        double diff = ch7[i] - ch8[i];
        te += k[i] * (h * diff);
    }
    const double te_max =
        std::max({ std::abs(te.r.x), std::abs(te.r.y), std::abs(te.r.z),
                   std::abs(te.v.x), std::abs(te.v.y), std::abs(te.v.z) });

    const double y_max =
        std::max({ std::abs(s.r.x), std::abs(s.r.y), std::abs(s.r.z),
                   std::abs(s.v.x), std::abs(s.v.y), std::abs(s.v.z) });
    const double te_allowed = std::max(y_max, 1.0) * tol;

    // Adaptive step size (1/8 exponent for 7th-order error)
    const double delta  = std::pow(te_allowed / (te_max + eps), 1.0 / 8.0);
    const double h_next = std::min(0.9 * delta * h, 4.0 * h);

    if (h_next < 16.0 * eps)
    {
        std::ostringstream ss;
        ss << "RKF78: next step fell below minimum at t=" << ti;
        throw AstroException(ss.str());
    }

    if (te_max > te_allowed)
    {
        // Step is rejected — return current state with reduced step
        return { s, et, TimeDelta(h_next), 0 };
    }

    // 8th-order solution
    PosState s_next = s;
    for (int i = 0; i < 13; ++i)
        s_next += k[i] * (h * ch8[i]);

    return { s_next, et + dt, TimeDelta(h_next), 0 };
}

std::vector<RKF78::Result> RKF78::doSteps(
    const ODE& ode, const PosState& s,
    const EphemerisTime& et0, const EphemerisTime& et1,
    const TimeDelta& dt)
{
    std::vector<Result> res;
    res.push_back({ s, et0, dt, 0 });

    while (res.back().et < et1)
    {
        res.push_back(doStep(ode, res.back().s, res.back().et, res.back().dt_next));
        if (res.back().et + res.back().dt_next > et1)
            res.back().dt_next = et1 - res.back().et;
    }

    return res;
}

} // namespace astro
