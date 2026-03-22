#include <iostream>

#include "State.h"
#include "Exceptions.h"
#include "SpiceCore.h"

#include <cspice/SpiceUsr.h>

namespace astro {

State State::transform(const ReferenceFrame& fromFr, const ReferenceFrame toFr, const EphemerisTime& et)
{
    if (fromFr.getType() == ReferenceFrameType::BodyFixedRotating &&
        toFr.getType()   == ReferenceFrameType::BodyFixedRotating)
    {
        throw AstroException("Either to-frame or from-frame must be inertial when transforming state");
    }

    if ((fromFr.isJ2000() || fromFr.getType() == ReferenceFrameType::BodyFixedNonRotating) &&
        (toFr.isJ2000()   || toFr.getType()   == ReferenceFrameType::BodyFixedNonRotating))
    {
        return *this; // No state transform needed
    }

    // Either to or from is inertial — get the body-fixed one
    Mat3 M;
    bool fromInertial;
    int body;
    if (toFr.isJ2000() || toFr.getType() == ReferenceFrameType::BodyFixedNonRotating)
    {
        // from body-fixed to inertial
        M            = fromFr.getRotationToJ2000(et);
        fromInertial = false;
        body         = fromFr.getCenterId();
    }
    else
    {
        // from inertial to body-fixed
        M            = glm::inverse(toFr.getRotationToJ2000(et));
        fromInertial = true;
        body         = toFr.getCenterId();
    }

    double tispm[6][6];
    {
        std::lock_guard<std::mutex> lock(Spice().mutex());
        tisbod_c("J2000", body, et.getETValue(), tispm);
    }

    if (!fromInertial)
    {
        double tispm_inv[6][6];
        invstm_c(tispm, tispm_inv);
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j)
                tispm[i][j] = tispm_inv[i][j];
    }

    // TODO: Retrieve angular velocities from the transformation
    std::cout << "WARNING - State::transform does not currently transform rotations" << std::endl;

    double state[6]    = { P.r.x, P.r.y, P.r.z, P.v.x, P.v.y, P.v.z };
    double tr_state[6] = {};
    mxvg_c(tispm, state, 6, 6, tr_state);

    State ret;
    ret.P.r = Vec3(tr_state[0], tr_state[1], tr_state[2]);
    ret.P.v = Vec3(tr_state[3], tr_state[4], tr_state[5]);
    return ret;
}


std::ostream& operator<<(std::ostream& os, const astro::PosState& s)
{
    os << "r: (" << s.r.x << ", " << s.r.y << ", " << s.r.z << ") [km]\n";
    os << "v: (" << s.v.x << ", " << s.v.y << ", " << s.v.z << ") [km/s]\n";
    return os;
}

PosState operator/(const PosState& p1, const PosState& p2)
{
    return PosState(p1.r / p2.r, p1.v / p2.v);
}

PosState abs(const PosState& p)
{
    return PosState(
        Vec3(std::abs(p.r.x), std::abs(p.r.y), std::abs(p.r.z)),
        Vec3(std::abs(p.v.x), std::abs(p.v.y), std::abs(p.v.z)));
}

std::ostream& operator<<(std::ostream& os, const astro::RotState& s)
{
    os << "q: (" << s.q.x << ", " << s.q.y << ", " << s.q.z << ", " << s.q.w << ") []\n";
    os << "w: (" << s.w.x << ", " << s.w.y << ", " << s.w.z << ") [RAD/s]\n";
    return os;
}

} // namespace astro
