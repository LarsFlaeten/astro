#ifndef _ASTRO_ORBIT_H_
#define _ASTRO_ORBIT_H_

// Contains base class for all orbits, as well as a simple orbit class for
// Kleplerian (2-body) orbits when the sattelite mass << orbit centre mass
// A lot of the algorithms here are taken from
// [1]  Orbital Mechanics for Engineering Students, 2nd Edition, Howard D. Curtis

// additionally taken from [2],[3]:
// https://downloads.rene-schwarz.com/download/M001-Keplerian_Orbit_Elements_to_Cartesian_State_Vectors.pdf
// https://downloads.rene-schwarz.com/download/M002-Cartesian_State_Vectors_to_Keplerian_Orbit_Elements.pdf 
#include "State.h"
#include "Time.h"
#include "Util.h"
#include <sstream>
#include <iostream>

namespace astro
{


// The base class off all orbits
class Orbit
{
public:
    Orbit();
    virtual ~Orbit();

    // The essential method of all orbits in astro:
    virtual State   getState(const EphemerisTime& et) = 0;

    // We need more stuff here, such as get frame etc...

};

// The six orbital elements of a Keplerian orbit
struct OrbitElements
{
    double h;   // Specific angular momentum
    double i;   // Inclination
    double omega; // RA of the ascending node
    double e;   // Eccentricity
    double w;   // Argument of perigee (originally this is small omega, but we use w...)
    double theta; // True anomaly at epoch
    double M0;  // Mean anomaly at epoch
    double rp;   // Periapsis distance
    double mu;   // The gravitaional parameter of the primary body
    EphemerisTime epoch;

    // Convert at state vector to orbit elements for the same frame of reference
    // Using method from [2]
    // mu is the gravitational parameter of the reference (attracting) body
    static OrbitElements fromStateVectorOE(const State& state, const EphemerisTime& epoch, double mu);
    
    	
    // Did some optimization, but only gained a few percent (and the code is not that clear anymore)
    // SKip it for now...
    //static OrbitElements fromstateVectorOEOpt(const State& state, double mu);

    // Using the Spice library
    static OrbitElements fromStateVectorSpice(const State& state,const EphemerisTime& epoch,  double mu);


    // Converts the current orbital elements to a state vector in the same frame of reference
    // Using [1]
    State   toStateVectorOE();


    // Using the spice method
    State   toStateVectorSpice();
};



template < class T >
std::ostream& operator << (std::ostream& os, const astro::OrbitElements& oe)
{
    double dpr = astro::DEGPERRAD;
    os << "Angular momentum:    " << oe.h       << " [km²/s]\n";
    os << "Inclination:         " << oe.i*dpr   << " [Deg]\n";
    os << "RA of the asc. node: " << oe.omega*dpr<< " [Deg]\n";
    os << "Eccentricity:        " << oe.e   << " [-]\n";
    os << "Argument of perigee: " << oe.w*dpr   << " [Deg]\n";   
    os << "True anomaly:        " << oe.theta*dpr<< " [Deg]\n";
    os << "Mean anomaly:        " << oe.M0*dpr  << " [Deg]\n";
    os << "Periapsis distance:  " << oe.rp      << " [km]\n"; 
    os << "Epoch:               " << oe.epoch.getETValue() << " [seconds]\n";
    os << "mu:                  " << oe.mu << " [km³/s²]\n";

    return os;
} 

void print(const astro::OrbitElements& oe)
{
    double dpr = astro::DEGPERRAD;
    std::cout << "Angular momentum:    " << oe.h       << " [km²/s]\n";
    std::cout << "Inclination:         " << oe.i*dpr   << " [Deg]\n";
    std::cout << "RA of the asc. node: " << oe.omega*dpr<< " [Deg]\n";
    std::cout << "Eccentricity:        " << oe.e   << " [-]\n";
    std::cout << "Argument of perigee: " << oe.w*dpr   << " [Deg]\n";
    std::cout << "True anomaly:        " << oe.theta*dpr<< " [Deg]\n";
    std::cout << "Mean anomaly:        " << oe.M0*dpr  << " [Deg]\n";
    std::cout << "Periapsis distance:  " << oe.rp      << " [km]\n";
    std::cout << "Epoch:               " << oe.epoch.getETValue() << " [seconds]\n";
    std::cout << "mu:                  " << oe.mu << " [km³/s²]\n";
    std::cout << std::endl; 


}


}



#endif
