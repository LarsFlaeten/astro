#include "Observer.h"
#include "SpiceCore.h"

namespace astro {

    Observer::Observer(int centerObj, ReferenceFrame refFrame, State state)
       : centerObj(centerObj), refFrame(refFrame), state(state)
    {
    } 
  
    State Observer::getState() const {
        return state;
    }
    
    ReferenceFrame Observer::getReferenceFrame() const {
        return refFrame;
    }

    int Observer::getCenterObject() const {
        return centerObj;
    }
    
    void    Observer::setState(const State& _state) {
        state = _state;
    }
            
    void    Observer::setCenterObject(int newCenter, bool recalcState, const EphemerisTime& et) {
        if(!refFrame.isJ2000())
            throw std::runtime_error("Change of center object can only be done in J2000 frame");

        if(recalcState) {
            // find the vector from the new center to the old:
            PosState c_state;
            Spice().getRelativeGeometricState(newCenter, centerObj, et, c_state, refFrame);
            state.P += c_state; 
            // TODO: Recalc orientation/rotState? Not if frame is J2000

        }
        centerObj = newCenter;
    }

    void    Observer::setReferenceFrame(ReferenceFrame newRF, bool recalcState, const EphemerisTime& et) {
        if(recalcState) {
            state = state.transform(refFrame, newRF, et);
        }                        

        refFrame = newRF;
    }


}
