#ifndef _ASTRO_UTIL_H
#define _ASTRO_UTIL_H

#include "Math.h"
#include <sstream>

namespace astro {
    // Pi and a few derived constants
    double pi();
    const double PI       = pi();
    const double TWOPI    = 2.0 * pi();
    const double PIHALF   = pi() * 0.5;

    // Degrees / radians conversion
    double degreesPerRadian();
    double radiansPerDegree();
    const double DEGPERRAD = degreesPerRadian();
    const double RADPERDEG = radiansPerDegree();

    // Convert to a position from Right Ascension, Declination
    void raDecToVec(double range, double ra, double dec, Vec3* res);

    // Convert from a position to Right Ascension, Declination
    void vecToRaDec(const Vec3& pos, double* range, double* ra, double* dec);

    // Normalizes any number to an arbitrary range by assuming the range wraps
    // around when going below min or above max
    double wrap(double value, double start, double end);

    // Transforms a vector by a quaternion sandwich product: q1 * pure(v) * q2
    // Usage (if Q denotes the orientation of a body):
    //   v_b = transform(Q_inv, v,   Q)
    //   v   = transform(Q,     v_b, Q_inv)
    Vec3 transform(const Quat& q1, const Vec3& v, const Quat& q2);

    // Pretty-printers
    std::ostream& operator<<(std::ostream& os, const Vec3& v);
    std::ostream& operator<<(std::ostream& os, const Quat& q);
    std::ostream& operator<<(std::ostream& os, const Mat3& m);

    // Orientation helpers
    Quat getProgradeOrientation(const Vec3& pos, const Vec3& velocity);
    Quat getRetrogradeOrientation(const Vec3& pos, const Vec3& velocity);

} // namespace astro

#endif
