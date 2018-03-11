#ifndef _ASTRO_STATE_H_
#define _ASTRO_STATE_H_

#include <sstream>

// Decided to use ork::math. Other users need to install ork for this to work.
// TODO: Fix dependency in CMakeLists so that this is clear for users
#ifndef ORK_API
#define ORK_API
#endif
#include <ork/math/vec3.h>


using ork::vec3d;

namespace astro
{


// The state of an object in a given reference frame
struct State {
    vec3d   r; // Position [km]
    vec3d   v; // Orbital Velocity [km/s]
};

// The derivative of a state vector
struct StateDot {
    vec3d   v; // Orbital velocity [km/s]
    vec3d   a; // Acceleration [km²/s]

};

std::ostream& operator << (std::ostream& os, const astro::State& s);

std::ostream& operator << (std::ostream& os, const astro::StateDot& dsdt);


}
#endif
