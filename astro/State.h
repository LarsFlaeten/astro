#ifndef _ASTRO_STATE_H_H
#define _ASTRO_STATE_H_H

#include <sstream>
#include <iostream>

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



template < class T >
std::ostream& operator << (std::ostream& os, const astro::State& s)
{
    os << "R: (" << s.r.x <<", " << s.r.y << ", " << s.r.z << ") [km]\n";
    os << "V: (" << s.v.x <<", " << s.v.y << ", " << s.v.z << ") [km/s]\n";
    return os;
}   

void print(const astro::State& s)
{
    std::cout << "R: (" << s.r.x <<", " << s.r.y << ", " << s.r.z << ") [km]\n";
    std::cout << "V: (" << s.v.x <<", " << s.v.y << ", " << s.v.z << ") [km/s]\n";
    std::cout << std::endl;
}

}
#endif
