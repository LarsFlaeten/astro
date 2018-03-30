#ifndef _ASTRO_ORBITELEMENTS_H_
#define _ASTRO_ORBITELEMENTS_H_

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
#include <utility>

namespace astro
{

const int KEPLER_MAX_ITERATIONS = 100;
const double KEPLER_TOLERANCE = 1.0E-10;


// Returns the asymptoe (angle between asymptote and the vector from perigee through locus)
// Undefined behavior for e < 1.0
double hyperbolicAsymptote(double eccentricty);

// Returns the hyperbolic excess velocity
// Undefined behavior for e < 1.0
double hyperbolicExcessVelocity(double mu, double a);





// The six orbital elements of a Keplerian orbit
struct OrbitElements
{
    // The five constant elements
    double h;   // Specific angular momentum
    double i;   // Inclination
    double omega; // RA of the ascending node
    double e;   // Eccentricity
    double w;   // Argument of perigee (originally this is small omega, but we use w...)
    
    // Time fix of this orbit
    double M0;  // Mean anomaly at epoch
    EphemerisTime epoch;
    
    // Auxillary variables
    double mu;   // The gravitaional parameter of the primary body
 
    // Derived quanitties:
    double rp;   // Periapsis distance
    double ap;   // Apoapsis distance
    double a;    // Semimajor axis 
    double T;    // Period (a negative number of e >= 1)
    double n;    // Mean motion

    // Computes the above 5 quantities when all thw others are given
    void computeDerivedQuantities();


	// Retrieves the state vector at the given time:
	// This method uses OE and not Spice, due to speed
	// (The two are otherwise equivalent)	
	State	toStateVector(const EphemerisTime& et);    	

    // Convert a state vector to orbit elements for the same frame of reference
	// at the given epoch
    // mu is the gravitational parameter of the reference (attracting) body
    static OrbitElements fromStateVector(
        const State& state,
        const EphemerisTime& epoch,
        double mu);



	// Below methods are mostly internal




    // Convert at state vector to orbit elements for the same frame of reference
    // Using method from [2]
    // mu is the gravitational parameter of the reference (attracting) body
    static OrbitElements fromStateVectorOE(
        const State& state,
        const EphemerisTime& epoch,
        double mu);

    // Using the Spice library
    static OrbitElements fromStateVectorSpice(
        const State& state,
        const EphemerisTime& epoch,
        double mu);


    // Given eccentricity e and true anomaly, computes the eccentric anomaly
    static double eccentricAnomalyFromTrueAnomaly(double trueAnomaly, double e);    

    // Given eccentricity e and mean anomaly, computes the eccentric anomaly
    static double eccentricAnomalyFromMeanAnomaly(double meanAnomaly, double e);
 

    // Given eccentricity e and true anomaly, computes the mean anomaly
    static double meanAnomalyFromTrueAnomaly(double trueAnomaly, double e);    

    // Given eccentricity e and mean anomaly, computes the true anomaly
    static double trueAnomalyFromMeanAnomaly(double meanAnomaly, double e);
    

    // Different Implementations of the solution of keplers equation
    // All functions take as input
    // double M - Mean anomaly
    // double e - Eccentricity
    // Return value: Ecentric Anomaly, number of iterations
    
    // Kepler1: Fixed point iteration
    static std::pair<double,int> Kepler1(double M, double e);
    // Kepler2: Newton Rapson
    static std::pair<double,int> Kepler2(double M, double e);
    // TODO? Implement even faster solver:
    // https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19720016564.pdf

    // Converts the current orbital elements to a state vector in the same
    // frame of reference. Using method from [1] ++
    State   toStateVectorOE(const EphemerisTime& et);


    // Converts the current orbital elements to a state vector in the same
    // frame of reference. Using the spice method
    // et: Time to resolve the state vector
    State   toStateVectorSpice(const EphemerisTime& et);

};




std::ostream& operator << (std::ostream& os, const astro::OrbitElements& oe);

void print(const astro::OrbitElements& oe);

}



#endif
