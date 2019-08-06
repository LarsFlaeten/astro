#ifndef _ASTRO_ABSERVER_H_
#define _ASTRO_ABSERVER_H_

#include "ReferenceFrame.h"
#include "State.h"
namespace astro {

    // An observer, which consist of a state relative to a spice object, and a reference
    // in which the state is given.
    // Observers are used to report positions of celestial objects relative to the observer, when
    // the observer state is not given by Spice, but calculated by Astro.
    // E.g planet positions relative to a space ship
    class Observer {
        public:
            Observer(int centerObj, ReferenceFrame refFrame, State state); 

            virtual State getState() const;
             
            virtual ReferenceFrame getReferenceFrame() const;

            virtual int getCenterObject() const;

            virtual void    setState(const State& state);

            // Sets a new center object for this observer frame
            // If recalcState is set to true, the state of the observer will be recalculated
            // for the new center object
            virtual void    setCenterObject(int newCenter, bool recalcState = false, const EphemerisTime& et = EphemerisTime(0));
            
            // Sets a new Reference Frame for this observers state
            // If recalcState is set to true, the state of the observer will be recalculated
            // for the new frame            
            virtual void    setReferenceFrame(ReferenceFrame rf, bool recalcState = false, const EphemerisTime& et = EphemerisTime(0));

        protected:
            int centerObj;
            ReferenceFrame refFrame;
            State   state;

    };



}


#endif

