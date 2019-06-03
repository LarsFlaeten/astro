#include "Orbit.h"
#include "OrbitElements.h"
#include "Util.h"
#include "SpiceCore.h"

#include <iostream>
#include <sstream>
#include <mutex>

#include <cspice/SpiceUsr.h>

#include <mork/math/vec3.h>
#include <mork/math/mat3.h>
#include <mork/math/mat4.h>

using mork::vec3d;
using mork::mat3d;
using mork::mat4d;

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

PosState SimpleOrbit::getState(const EphemerisTime& et)
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
