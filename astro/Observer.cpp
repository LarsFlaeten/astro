#include "Observer.h"
#include "SpiceCore.h"

namespace astro {

    Observer::Observer(int centerObj, ReferenceFrame refFrame, PosState posState)
       : centerObj(centerObj), refFrame(refFrame), posState(posState)
    {
    } 
  
    PosState Observer::getState() const {
        return posState;
    }
    
    ReferenceFrame Observer::getReferenceFrame() const {
        return refFrame;
    }

    int Observer::getCenterObject() const {
        return centerObj;
    }
    
    void    Observer::setState(const PosState& state) {
        posState = state;
    }
            
    void    Observer::setCenterObject(int newCenter, bool recalcState, const EphemerisTime& et) {
        if(recalcState) {
            // find the vector from the new center to the old:
            PosState c_state;
            Spice().getRelativeGeometricState(newCenter, centerObj, et, c_state, refFrame);
            posState += c_state; 


        }
        centerObj = newCenter;
    }

    void    Observer::setReferenceFrame(ReferenceFrame newRF, bool recalcState, const EphemerisTime& et) {
        if(recalcState) {
            PosState newState = posState.transform(refFrame, newRF, et);
            posState = newState;
        }                        

        refFrame = newRF;
    }


}
