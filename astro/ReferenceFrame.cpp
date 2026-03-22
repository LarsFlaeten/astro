#include "ReferenceFrame.h"
#include "SpiceCore.h"

#include <cspice/SpiceUsr.h>
#include <sstream>

namespace astro {

ReferenceFrame::ReferenceFrame()
{}

ReferenceFrame::~ReferenceFrame()
{}

Mat3 ReferenceFrame::getRotationToJ2000(const EphemerisTime& et) const
{
    if (type == ReferenceFrameType::Inertial || type == ReferenceFrameType::BodyFixedNonRotating)
        return Mat3(1.0);

    // Use Spice to get the body-fixed → inertial rotation.
    // tipbod_c returns tipm such that v_body = tipm * v_inertial.
    // We want the inverse (inertial → body), so we return its inverse.
    double tipm[3][3];
    {
        std::lock_guard<std::mutex> lock(Spice().mutex());
        tipbod_c("J2000", centerId, et.getETValue(), tipm);
    }
    Spice().checkError();

    // SPICE tipm is row-major; construct GLM column-major matrix.
    Mat3 tip;
    for (int col = 0; col < 3; ++col)
        for (int row = 0; row < 3; ++row)
            tip[col][row] = tipm[row][col];

    return glm::inverse(tip);
}

ReferenceFrameType ReferenceFrame::getType() const
{
    return type;
}

int ReferenceFrame::getId() const
{
    return spiceId;
}

int ReferenceFrame::getCenterId() const
{
    return centerId;
}

std::string ReferenceFrame::getName() const
{
    return spiceName;
}

bool ReferenceFrame::operator==(const ReferenceFrame& other) const
{
    return centerId == other.centerId &&
           type     == other.type     &&
           spiceId  == other.spiceId;
}

ReferenceFrame ReferenceFrame::createJ2000()
{
    // Manually create this frame: Spice returns the wrong frame for ref id 0
    // (it returns the old MARSIAU frame).
    ReferenceFrame ref;
    ref.type      = Inertial;
    ref.spiceName = "J2000";
    ref.spiceId   = 1;
    ref.centerId  = 0; // SSB
    return ref;
}

ReferenceFrame ReferenceFrame::fromString(const std::string& name)
{
    if (name == "J2000")
        return createJ2000();
    throw std::runtime_error("Frame name '" + name + "' is not implemented in fromString()");
}

bool ReferenceFrame::isJ2000() const
{
    return type == Inertial && spiceId == 1 && centerId == 0;
}

ReferenceFrame ReferenceFrame::createBodyFixedSpice(int bodyId)
{
    if (bodyId == 0)
        return ReferenceFrame::createJ2000();

    ReferenceFrame ref;
    ref.type     = (bodyId < 10) ? BodyFixedNonRotating : BodyFixedRotating;
    ref.centerId = bodyId;

    {
        std::lock_guard<std::mutex> lock(Spice().mutex());
        const int lenout = 32;
        char  frameName[lenout];
        int   frameId;
        int   found;
        cidfrm_c(bodyId, lenout, &frameId, frameName, &found);

        if (!found)
        {
            std::ostringstream ss;
            ss << "Cannot create body-fixed frame for id " << bodyId
               << ": object id not found in loaded kernels";
            throw std::runtime_error(ss.str());
        }
        ref.spiceName = frameName;
        ref.spiceId   = frameId;
    }
    Spice().checkError();

    return ref;
}

} // namespace astro
