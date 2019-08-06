#include "ReferenceFrame.h"
#include "SpiceCore.h"
#include <cspice/SpiceUsr.h>


namespace astro {

ReferenceFrame::ReferenceFrame() {
}

ReferenceFrame::~ReferenceFrame() {
}


mork::mat3d  ReferenceFrame::getRotationToJ2000(const EphemerisTime& et) const {

    // If its the J2000, return identity:
    if(this->getType() == ReferenceFrameType::Inertial || this->getType() == ReferenceFrameType::BodyFixedNonRotating)
        return mork::mat3d::IDENTITY;

    mork::mat3d tip;
    // If not, we use spice to get the relative rotation:
    {
        std::lock_guard<std::mutex> lock(Spice().mutex());
        // We use tipbod:
        // Return a 3x3 matrix that transforms positions in inertial 
        // coordinates to positions in body-equator-and-prime-meridian 
        // coordinates. 
        // ie v_body = M * v_inertial;         
        // (We want the inverse of this)
        double tipm[3][3];
        tipbod_c("J2000", centerId, et.getETValue(), tipm);
        tip = mork::mat3d(tipm);

        //for(int i = 0; i < 3; ++i)
        //    std::cout << "[ " << tipm[i][0] << ", " << tipm[i][1] << ", " << tipm[i][2] << " ]\n";
        //std::cout << std::endl;

    }
    Spice().checkError();
 


    return tip.inverse();


}

ReferenceFrameType ReferenceFrame::getType() const {
    return type;
}

int ReferenceFrame::getId() const {
    return spiceId;
}

int ReferenceFrame::getCenterId() const {
    return centerId;
}

std::string ReferenceFrame::getName() const {
    return spiceName;
}

bool ReferenceFrame::operator==(const ReferenceFrame& other) const {
    if(this->centerId == other.centerId &&
            this->type == other.type &&
            this->spiceId == other.spiceId)
        return true;
    else
        return false;

}


ReferenceFrame  ReferenceFrame::createJ2000() {
    
    // We manually create this frame with name and id,
    // as spice will return wrong frame when using ref id 0 (it return the old
    // MARSIAU frame
    ReferenceFrame ref;
    ref.type = Inertial;
    ref.spiceName = "J2000";
    ref.spiceId = 1;
    ref.centerId = 0; // SSB
    return ref;
}

ReferenceFrame  ReferenceFrame::fromString(const std::string& name) {
    
    // We manually create this frame with name and id,
    // as spice will return wrong frame when using ref id 0 (it return the old
    // MARSIAU frame
    ReferenceFrame ref;

    if(name.compare("J2000") == 0) {
        ref.type = Inertial;
        ref.spiceName = "J2000";
        ref.spiceId = 1;
        ref.centerId = 0; // SSB
    } else
        throw std::runtime_error("Frame name " + name + " is not implemented in fromString() factory");

    return ref;
}


bool    ReferenceFrame::isJ2000() const {
    return (type == Inertial && spiceId == 1 && centerId == 0);

}

ReferenceFrame ReferenceFrame::createBodyFixedSpice(int bodyId) {
    if(bodyId == 0)
        return ReferenceFrame::createJ2000();
    
    ReferenceFrame ref;
    if(bodyId < 10)
        ref.type = BodyFixedNonRotating;
    else
        ref.type = BodyFixedRotating;

    ref.centerId = bodyId;
    // Try to get the frame name and id from Spice. If appropriate kernels are not loaded,
    // this will throw an exception
    {
        std::lock_guard<std::mutex> lock(Spice().mutex());
        
        int lenout = 32;
        char frameName[lenout];
        int frameId;
        int found;
        cidfrm_c(bodyId, lenout, &frameId, frameName, &found);

        if(!found) {
            std::stringstream ss;
            ss << "Cannot create body fixed frame for id " << bodyId << ", object id not found in loaded kernels";
            
            throw std::runtime_error(ss.str());
        
        }
        ref.spiceName = frameName;
        ref.spiceId = frameId;    


    }
    Spice().checkError();
     

    return ref;

}

}
