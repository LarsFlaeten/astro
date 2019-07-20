#ifndef _ASTRO_REFERENCE_FRAME_H_
#define _ASTRO_REFERENCE_FRAME_H_

#include <mork/math/mat3.h>

#include "Time.h"

namespace astro {

enum ReferenceFrameType {
    Inertial,               // Allways J2000 axes
    BodyFixedNonRotating, // Body fixed, but same axis as inertial (J2000)
    BodyFixedRotating,  // Rotating with the body

};


class ReferenceFrame {

public:
    ReferenceFrame();

    virtual ~ReferenceFrame();

    // Returns the orientation/rotation of the frame in J2000 inertial frame
    // This matrix can be used to get the orientation of a body in inertial frame,
    // or to transform body fixed coordinates to inertial:
    // r_j2000 = M * r_body. To transform from inertial to body fixed,
    // use M.inverse()
    virtual mork::mat3d  getRotation(const EphemerisTime& et) const;

    virtual ReferenceFrameType getType() const;

    // returns the Id of the center object of this frame
    virtual int getCenterId() const;
    
    // Returns the internal (spice) ID of this frame
    virtual int getId() const;

   
    virtual std::string getName() const;

    virtual bool operator==(const ReferenceFrame& other) const;

    // Creates a J2000 inertial frame. This is the same as the ICRF inertial frame.
    static ReferenceFrame createJ2000();

    // Creates a body fixed frame with bodyId as center object
    // For planets and sattelites, this will be a rotating frame
    // For barycenters, it will be a non-rotating frame
    static ReferenceFrame createBodyFixedSpice(int bodyId);


    // Returns true of this is the J2000 frame
    virtual bool isJ2000() const;

protected:
    ReferenceFrameType type;

    int spiceId;
    std::string spiceName;
    int centerId;
};


}


#endif
