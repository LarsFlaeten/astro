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
    Orbit();
    virtual ~Orbit();

    // The essential method of all orbits in astro:
    virtual PosState   getState(const EphemerisTime& et) = 0;

    // We need more stuff here, such as get frame etc...


};

class SimpleOrbit : public Orbit
{
public:
    SimpleOrbit(const OrbitElements& oe);
    virtual ~ SimpleOrbit();

    virtual PosState getState(const EphemerisTime& et);

    // returns the period of the orbit. If this is not a periodic orbit
    // (parabolic, hyperbolic) the return value is negative
    // Return value: Period [seconds]
    // TODO: Evaluate to move this up to Base class
    // TODO: Evaluate to throw exception for non-repeating/open orbits
    double  getPeriod() const;

    // Returns true if this is an elliptic or circular orbit, false otherwise
    bool    isPeriodic() const; 

    // Returns a reference to te underøying orbital elements of this Orbit
    const OrbitElements& getOrbitElements() const;

protected:
    SimpleOrbit();


    OrbitElements oe;

};

}

#endif
