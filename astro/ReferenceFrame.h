#ifndef _ASTRO_REFERENCE_FRAME_H_
#define _ASTRO_REFERENCE_FRAME_H_

#include "Math.h"
#include "Time.h"

namespace astro {

enum ReferenceFrameType
{
    Inertial,             // Always J2000 axes
    BodyFixedNonRotating, // Body fixed, same axes as inertial (J2000)
    BodyFixedRotating,    // Rotating with the body
};


class ReferenceFrame
{
public:
    ReferenceFrame();
    virtual ~ReferenceFrame();

    // Returns the orientation/rotation of the frame in J2000 inertial frame.
    // r_j2000 = M * r_body. To transform from inertial to body fixed, use
    // glm::inverse(M).
    virtual Mat3 getRotationToJ2000(const EphemerisTime& et) const;

    virtual ReferenceFrameType getType() const;

    // Returns the Id of the center object of this frame
    virtual int getCenterId() const;

    // Returns the internal (Spice) ID of this frame
    virtual int getId() const;

    virtual std::string getName() const;

    virtual bool operator==(const ReferenceFrame& other) const;

    // Creates a J2000 inertial frame (same as ICRF)
    static ReferenceFrame createJ2000();

    // Creates a body-fixed frame with bodyId as center object.
    // For planets and satellites this is a rotating frame;
    // for barycenters it is non-rotating.
    static ReferenceFrame createBodyFixedSpice(int bodyId);

    // Creates a frame based on its name. Currently only "J2000" is implemented.
    static ReferenceFrame fromString(const std::string& name);

    virtual bool isJ2000() const;

protected:
    ReferenceFrameType type;
    int spiceId;
    std::string spiceName;
    int centerId;
};

} // namespace astro

#endif
