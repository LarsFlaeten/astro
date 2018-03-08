#ifndef _ASTRO_ORBIT_H_
#define _ASTRO_ORBIT_H_

// Contains base class for all orbits, as well as a simple orbit class for
// Kleplerian (2-body) orbits when the sattelite mass << orbit centre mass

#include "State.h"
#include "Time.h"
#include "Util.h"
#include "OrbitElements.h"
#include <sstream>
#include <iostream>
#include <utility>

namespace astro
{

// The base class off all orbits
class Orbit
{
public:
    Orbit(const OrbitElements& oe);
    virtual ~Orbit();

    // The essential method of all orbits in astro:
    virtual State   getState(const EphemerisTime& et) = 0;

    // We need more stuff here, such as get frame etc...

protected:
    Orbit();

};


}

#endif
