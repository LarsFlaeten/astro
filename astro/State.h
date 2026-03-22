#ifndef _ASTRO_STATE_H_
#define _ASTRO_STATE_H_

#include <sstream>

#include "Math.h"
#include "ReferenceFrame.h"
#include "Time.h"

namespace astro {

// The translative state of an object in a given reference frame
class PosState
{
public:
    Vec3 r; // Position [km]
    Vec3 v; // Orbital Velocity [km/s]

    PosState()
        : r(0.0), v(0.0)
    {}

    PosState(const Vec3& _r, const Vec3& _v)
        : r(_r), v(_v)
    {}

    explicit PosState(double val)
        : r(val), v(val)
    {}

    PosState& operator+=(const PosState& o)
    {
        r += o.r;
        v += o.v;
        return *this;
    }

    PosState& operator-=(const PosState& o)
    {
        r -= o.r;
        v -= o.v;
        return *this;
    }

    PosState& operator*=(double a)
    {
        r *= a;
        v *= a;
        return *this;
    }

    friend PosState operator+(PosState lhs, const PosState& rhs) { return lhs += rhs; }
    friend PosState operator-(PosState lhs, const PosState& rhs) { return lhs -= rhs; }
    friend PosState operator*(PosState lhs, double a)            { return lhs *= a; }
    friend PosState operator*(double a, PosState rhs)            { return rhs *= a; }
};

// Only required for steppers with error control
PosState operator/(const PosState& p1, const PosState& p2);
PosState abs(const PosState& p);

std::ostream& operator<<(std::ostream& os, const astro::PosState& s);


// The rotation state of an object in a given reference coordinate system.
// Usually given in the global/inertial frame.
class RotState
{
public:
    Quat q; // Rotation quaternion in the given frame
    Vec3 w; // Angular velocities in the given frame [radians/s]

    RotState()
        : q(1.0, 0.0, 0.0, 0.0), w(0.0)
    {}

    RotState(const Quat& _q, const Vec3& _w)
        : q(_q), w(_w)
    {}
};

std::ostream& operator<<(std::ostream& os, const astro::RotState& s);


// A state of an object. Fundamentally the state consists of two parts
// which are uncoupled: translation (PosState) and rotation (RotState).
// The state does not know anything about its frame of reference.
class State
{
public:
    State() : P(), R() {}
    State(const PosState& p, const RotState& r) : P(p), R(r) {}

    PosState P;
    RotState R;

    // Transforms this state between reference frames at the given ET
    State transform(const ReferenceFrame& fromFr, const ReferenceFrame toFr, const EphemerisTime& et);
};

} // namespace astro

#endif
