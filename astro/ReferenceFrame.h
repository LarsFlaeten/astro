#ifndef _ASTRO_REFERENCE_FRAME_H_
#define _ASTRO_REFERENCE_FRAME_H_

#ifndef ORK_API
#define ORK_API
#endif
#include <ork/math/mat3.h>

#include "Time.h"

namespace astro {

enum ReferenceFrameType {
    Inertial,
    BodyFixed

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
    virtual ork::mat3d  getRotation(const EphemerisTime& et) const;

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
    static ReferenceFrame createBodyFixedSpice(int bodyId);



protected:
    ReferenceFrameType type;

    int spiceId;
    std::string spiceName;
    int centerId;
};


}


#endif
