#include "OrbitElements.h"
#include "Util.h"
#include "SpiceCore.h"

#include <iostream>
#include <sstream>
#include <mutex>

#include <cspice/SpiceUsr.h>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/transform.hpp>

namespace astro {

double hyperbolicAsymptote(double eccentricity)
{
    return PI - acos(-1.0 / eccentricity);
}

double hyperbolicExcessVelocity(double mu, double a)
{
    return std::sqrt(-mu / a);
}

OrbitElements OrbitElements::fromStateVector(const PosState& state, const EphemerisTime& epoch, double mu)
{
    return OrbitElements::fromStateVectorOE(state, epoch, mu);
}


OrbitElements OrbitElements::fromStateVectorOE(const PosState& state, const EphemerisTime& epoch, double mu)
{
    if (mu < 0)
        throw AstroException("Non-positive mass given for primary body");

    OrbitElements oe;

    // 1. Distance
    double r = glm::length(state.r);
    if (r < 1.0E-5)
        throw AstroException("Radius of state vector near zero — degenerate case");

    // 2. Scalar speed
    double v = glm::length(state.v);
    if (v < 1.0E-5)
        throw AstroException("Velocity near zero — degenerate case");

    // 3. Radial velocity: v_r = (R · V) / r
    double v_r = glm::dot(state.r, state.v / r);

    // 4. Specific angular momentum: H = R × V
    Vec3 H = glm::cross(state.r, state.v);

    // 5. Magnitude of angular momentum
    oe.h = glm::length(H);
    if (oe.h < 1.0E-5)
        throw AstroException("Angular momentum near zero — degenerate case");

    // 6. Inclination
    oe.i = std::acos(H.z / oe.h);

    // 7. Node line: N = K × H, K is Z-unit vector
    Vec3 N(-H.y, H.x, 0.0);

    // 8. Magnitude of node line
    double n = glm::length(N);

    // 9. RA of the ascending node
    oe.omega = (n > 0.0) ? std::acos(N.x / n) : 0.0;
    if (N.y < 0.0)
        oe.omega = 2.0 * astro::PI - oe.omega;

    // 10. Eccentricity vector
    Vec3 R_prime = state.r * (v * v - mu / r);
    Vec3 V_prime = state.v * (r * v_r);
    Vec3 E       = (R_prime - V_prime) * (1.0 / mu);

    // 11. Eccentricity
    oe.e = glm::length(E);

    // 12. Argument of perigee
    Vec3 E_norm = glm::normalize(E);
    if (n > 0.0)
    {
        Vec3 N_norm = glm::normalize(N);
        oe.w = std::acos(glm::dot(N_norm, E_norm));
    }
    else
    {
        oe.w = 0.0;
    }
    if (E.z < 0.0)
        oe.w = 2.0 * astro::PI - oe.w;

    // 13. True anomaly → mean anomaly
    Vec3   R_norm = glm::normalize(state.r);
    double theta  = std::acos(glm::dot(E_norm, R_norm));
    if (v_r < 0.0)
        theta = 2.0 * astro::PI - theta;

    oe.M0   = meanAnomalyFromTrueAnomaly(theta, oe.e);
    oe.mu   = mu;
    oe.epoch = epoch;
    oe.computeDerivedQuantities();

    return oe;
}

void OrbitElements::computeDerivedQuantities()
{
    rp = h * h / mu * (1.0 / (1.0 + e));
    ap = h * h / mu * (1.0 / (1.0 - e));
    a  = 0.5 * (rp + ap);
    T  = (e < 1.0) ? astro::TWOPI / std::sqrt(mu) * std::pow(a, 1.5) : -1.0;
    n  = std::sqrt(mu / std::pow(std::abs(a), 3.0));
}

double OrbitElements::eccentricAnomalyFromTrueAnomaly(double trueAnomaly, double e)
{
    if (e < 0.0)
        throw astro::AstroException("ERROR, negative eccentricity is not allowed");

    if (e < 1.0)
    {
        return 2.0 * std::atan(std::tan(trueAnomaly / 2.0) / std::sqrt((1 + e) / (1 - e)));
    }
    else
    {
        // Hyperbolic anomaly
        // http://control.asu.edu/Classes/MAE462/462Lecture05.pdf
        return 2.0 * std::atanh(std::sqrt((e - 1) / (e + 1)) * std::tan(trueAnomaly / 2.0));
    }
}

double OrbitElements::meanAnomalyFromTrueAnomaly(double trueAnomaly, double e)
{
    if (e < 0.0)
        throw astro::AstroException("ERROR, negative eccentricity is not allowed");

    if (e < 1.0)
    {
        double E = eccentricAnomalyFromTrueAnomaly(trueAnomaly, e);
        return E - e * std::sin(E);
    }
    else
    {
        double H = eccentricAnomalyFromTrueAnomaly(trueAnomaly, e);
        return e * std::sinh(H) - H;
    }
}

double OrbitElements::eccentricAnomalyFromMeanAnomaly(double M, double e)
{
    if (e < 0.0)
        throw astro::AstroException("ERROR, negative eccentricity is not allowed");
    return Kepler2(M, e).first;
}

double OrbitElements::trueAnomalyFromMeanAnomaly(double M, double e)
{
    if (e < 0.0)
        throw astro::AstroException("ERROR, negative eccentricity is not allowed");

    if (e < 0.9998)
    {
        double E     = eccentricAnomalyFromMeanAnomaly(M, e);
        double cosE2 = std::cos(0.5 * E);
        double sinE2 = std::sin(0.5 * E);
        return 2.0 * std::atan2(std::sqrt(1.0 + e) * sinE2, std::sqrt(1.0 - e) * cosE2);
    }
    else
    {
        double H = eccentricAnomalyFromMeanAnomaly(M, e);
        return 2.0 * std::atan(std::sqrt((e + 1.0) / (e - 1.0)) * std::tanh(0.5 * H));
    }
}


std::pair<double, int> OrbitElements::Kepler1(double M, double e)
{
    if (e < 0.0)
        throw astro::AstroException("ERROR, negative eccentricity is not allowed");
    if (e >= 1.0)
        throw astro::AstroException("ERROR, Kepler1 is for e < 1.0");

    double Mw = astro::wrap(M, 0.0, astro::TWOPI);
    double E0 = Mw;
    double En;
    int it = 0;

    while (true)
    {
        ++it;
        En = Mw + e * std::sin(E0);
        if (std::abs(En - E0) < 1.0E-10)
            break;
        if (it >= KEPLER_MAX_ITERATIONS)
            throw astro::AstroException("Kepler1 did not converge within max iterations");
        E0 = En;
    }
    return { En, it };
}

std::pair<double, int> OrbitElements::Kepler2(double M, double e)
{
    if (e < 0.0)
        throw astro::AstroException("ERROR, negative eccentricity is not allowed");

    if (e < 0.9998)
    {
        double Mw = astro::wrap(M, 0.0, astro::TWOPI);
        double E0 = Mw;
        double En;
        int it = 0;

        while (true)
        {
            ++it;
            En = E0 - (E0 - e * std::sin(E0) - Mw) / (1.0 - e * std::cos(E0));

            if (std::abs(En - E0) < KEPLER_TOLERANCE)
                break;
            if (it >= KEPLER_MAX_ITERATIONS)
            {
                std::ostringstream oss;
                oss << "Kepler2 did not converge within max iterations ("
                    << astro::KEPLER_MAX_ITERATIONS << "), e=" << e << ", M=" << Mw;
                throw astro::AstroException(oss.str());
            }
            E0 = En;
        }
        return { En, it };
    }
    else if (e > 1.0)
    {
        double H0 = M;
        double Hn;
        int it = 0;

        while (true)
        {
            ++it;
            Hn = H0 + (M - e * std::sinh(H0) + H0) / (e * std::cosh(H0) - 1.0);

            if (std::abs(Hn - H0) < KEPLER_TOLERANCE)
                break;
            if (it >= KEPLER_MAX_ITERATIONS)
            {
                std::ostringstream oss;
                oss << "Kepler2 did not converge within max iterations ("
                    << astro::KEPLER_MAX_ITERATIONS << "), e=" << e << ", M=" << M;
                throw astro::AstroException(oss.str());
            }
            H0 = Hn;
        }
        return { Hn, it };
    }
    else
    {
        throw astro::AstroException("ERROR, Kepler2 is for e < 0.9998 or e > 1.0");
    }
}

PosState OrbitElements::toStateVector(const EphemerisTime& et)
{
    return toStateVectorOE(et);
}

PosState OrbitElements::toStateVectorOE(const EphemerisTime& et)
{
    double t0    = epoch.getETValue();
    double t     = et.getETValue();
    double M     = M0 + n * (t - t0);
    double theta = trueAnomalyFromMeanAnomaly(M, e);

    double costheta = std::cos(theta);
    double sintheta = std::sin(theta);

    Vec3 Rxp(costheta, sintheta, 0.0);
    Rxp *= (h * h) / mu * (1.0 / (1.0 + e * costheta));

    Vec3 Vxp(-sintheta, e + costheta, 0.0);
    Vxp *= mu / h;

    // Rotation matrices from perifocal to geocentric equatorial frame.
    // Using glm::rotate on a mat4 identity and casting to mat3.
    // Note: the transpose is taken below to match the reference convention — see [1].
    Mat3 R3_w     = Mat3(glm::rotate(glm::dmat4(1.0), -w,     glm::dvec3(0, 0, 1)));
    Mat3 R1_i     = Mat3(glm::rotate(glm::dmat4(1.0), -i,     glm::dvec3(1, 0, 0)));
    Mat3 R3_Omega = Mat3(glm::rotate(glm::dmat4(1.0), -omega, glm::dvec3(0, 0, 1)));

    Mat3 Q_X_to_xp = R3_w * R1_i * R3_Omega;
    Mat3 Q_xp_to_X = glm::transpose(Q_X_to_xp);

    PosState state;
    state.r = Q_xp_to_X * Rxp;
    state.v = Q_xp_to_X * Vxp;
    return state;
}

OrbitElements OrbitElements::fromStateVectorSpice(const PosState& state, const EphemerisTime& epoch, double mu)
{
    astro::Spice(); // ensure Spice is initialised

    double et   = epoch.getETValue();
    double elts[8];

    {
        std::lock_guard<std::mutex> lock(astro::Spice().mutex());
        oscelt_c(const_cast<double*>(&(state.r.x)), et, mu, elts);
    }
    astro::Spice().checkError();

    OrbitElements oe;
    oe.rp    = elts[0];
    oe.e     = elts[1];
    oe.i     = elts[2];
    oe.omega = elts[3];
    oe.w     = elts[4];
    oe.h     = glm::length(glm::cross(state.r, state.v));
    oe.M0    = elts[5];
    oe.mu    = mu;
    oe.epoch = epoch;
    oe.computeDerivedQuantities();

    return oe;
}

PosState OrbitElements::toStateVectorSpice(const EphemerisTime& et)
{
    double elts[8] = { rp, e, i, omega, w, M0, epoch.getETValue(), mu };
    double st[6];
    {
        std::lock_guard<std::mutex> lock(astro::Spice().mutex());
        conics_c(elts, et.getETValue(), st);
    }
    astro::Spice().checkError();

    PosState state;
    state.r = Vec3(st[0], st[1], st[2]);
    state.v = Vec3(st[3], st[4], st[5]);
    return state;
}


std::ostream& operator<<(std::ostream& os, const astro::OrbitElements& oe)
{
    double dpr = astro::DEGPERRAD;
    os << "Angular momentum:    " << oe.h            << " [km²/s]\n";
    os << "Inclination:         " << oe.i * dpr      << " [Deg]\n";
    os << "RA of the asc. node: " << oe.omega * dpr  << " [Deg]\n";
    os << "Eccentricity:        " << oe.e             << " [-]\n";
    os << "Argument of perigee: " << oe.w * dpr       << " [Deg]\n";
    os << "Mean anomaly @epoch: " << oe.M0 * dpr      << " [Deg]\n";
    os << "Epoch:               " << oe.epoch.getETValue() << " [seconds]\n";
    os << "mu:                  " << oe.mu             << " [km³/s²]\n";
    os << "Periapsis distance:  " << oe.rp             << " [km]\n";
    os << "Apoapsis distance:   " << oe.ap             << " [km]\n";
    os << "Semimajor axis:      " << oe.a              << " [km]\n";
    os << "Period:              " << oe.T              << " [s]\n";
    os << "Mean motion:         " << oe.n              << " [rad/s]\n";
    return os;
}

void print(const astro::OrbitElements& oe)
{
    std::cout << oe;
    std::cout << std::endl;
}

} // namespace astro
