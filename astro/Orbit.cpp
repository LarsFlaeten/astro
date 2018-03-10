#include "Orbit.h"
#include "OrbitElements.h"
#include "Util.h"
#include "SpiceCore.h"

#include <iostream>
#include <sstream>
#include <mutex>

#include <cspice/SpiceUsr.h>

#ifndef ORK_API
#define ORK_API
#endif
#include <ork/math/vec3.h>
#include <ork/math/mat3.h>
#include <ork/math/mat4.h>

using ork::vec3d;
using ork::mat3d;
using ork::mat4d;

namespace astro
{


Orbit::Orbit()
{

}

Orbit::~Orbit()
{

}



SimpleOrbit::SimpleOrbit()
{

}

SimpleOrbit::~SimpleOrbit()
{

}

SimpleOrbit::SimpleOrbit(const OrbitElements& orbitElements)
    : Orbit(), oe(orbitElements)
{
    // Make sure all derived elements are calculated
    oe.computeDerivedQuantities();
}

State SimpleOrbit::getState(const EphemerisTime& et)
{
    return oe.toStateVector(et);
}

double SimpleOrbit::getPeriod() const
{
    return oe.T;
}

bool SimpleOrbit::isPeriodic() const
{
    return oe.e < 1.0 ? true : false;
}

const OrbitElements& SimpleOrbit::getOrbitElements() const
{
    return oe;
}

}
