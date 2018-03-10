#include <iostream>

#include <cxxopts.hpp>


#include <astro/SpiceCore.h>
#include <astro/Time.h>
#include <astro/State.h>
#include <astro/Orbit.h>


int main(int argc, char **argv)
{
    int example;
    cxxopts::Options options(argv[0], "Astrodynamics Library\n(c) 2018 Lars Flæten");
    
    try
    {
        options.add_options()
            ("h,help", "Print help")
            ("e,example", "Choose example", cxxopts::value<int>());


        auto result = options.parse(argc, argv);
        if(result.count("help"))
        {
            std::cout << options.help({""}) << std::endl;
            return 0;
        }

        if(result.count("example"))
        {
            example = result["example"].as<int>();
        } else
            std::cout << "Please give an example to run with -e <example>" << std::endl;
    }
    catch(cxxopts::OptionException){
        std::cout << options.help({""}) << std::endl;
        return -2;
        
    }
   // **********************************************************
    // Example 1 - TIME
    // This example shows normal time access and manipulation
    // **********************************************************
    if(example == 1)
    {
        // Load the lepseconds kernel
        astro::Spice().loadKernel("../data/spice/lsk/naif0012.tls");

        // The basic time keeping entity used in astro is Ephemeris time,
        // Same as Spice. This value is represented as seconds past the
        // J2000 epoch in the time system known as Barycentric Dynamical Time
        // (TDB).
        //
        // Some examples of initialization
        astro::EphemerisTime et1; // et = 0 represents the J2000 epoch
    
        std::cout << "Date/Time in ISO: " << et1.toISOUTCString() << ", in Julian Day(UTC): " << et1.toJDUTCString() << std::endl;
    
        et1 = astro::EphemerisTime::fromString("2018 February 22, 20:04:00 UTC"); // The time this line was written
        std::cout << "Date/Time in ISO: " << et1.toISOUTCString() << ", in Julian Day(UTC): " << et1.toJDUTCString() << std::endl;
    
        et1 = astro::EphemerisTime::fromString("2017-12-24T17:00:00.12"); // ISO format
        std::cout << "Date/Time in ISO: " << et1.toISOUTCString() << ", in Julian Day(UTC): " << et1.toJDUTCString() << std::endl;
    
        et1 = astro::EphemerisTime::fromString("2451515.2981 JD"); // From Julian date
        std::cout << "Date/Time in ISO: " << et1.toISOUTCString() << ", in Julian Day(UTC): " << et1.toJDUTCString() << std::endl;
 
        // Init with J2000, shall give 2000 JAN 01 12:00:00
        et1 = astro::EphemerisTime::fromJDUTC(2451545.0);
        std::cout << "Date/Time in ISO: " << et1.toISOUTCString() << ", in Julian Day(UTC): " << et1.toJDUTCString() << std::endl;
    }

    
    // **************************************************************************
    // Example 2 - Simple Orbit & Orbit Elements
    // This example shows how to establish a keplerian orbit from a state vector
    // and dumps the orbital elements
    // **************************************************************************
    if(example == 2)
    {
        astro::State   state;
        state.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
        state.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]

        astro::EphemerisTime et(0); // et = 0 represents the J2000 epoch
        double mu_earth = 398600.0;
    
        // Convert to Keplerian orbital elements for this epoch
        astro::OrbitElements oe = astro::OrbitElements::fromStateVector(state, et, mu_earth);
 
        // Generate a simple orbit from these elements:
        astro::SimpleOrbit orbit1(oe);    
    
        std::cout << orbit1.getOrbitElements() << std::endl;
         
        
    }

    // **************************************************************************
    // Example 3 - Simple Orbit Plot
    // This example shows how to establish a keplerian orbit from a state vector
    // and plots the position as a function of time
    // **************************************************************************
    if(example == 3)
    {
        astro::State   state;
        state.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
        state.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]

        astro::EphemerisTime et(0); // et = 0 represents the J2000 epoch
        double mu_earth = 398600.0;
    
        // Convert to Keplerian orbital elements for this epoch
        astro::OrbitElements oe = astro::OrbitElements::fromStateVector(state, et, mu_earth);
 
        // Generate a simple orbit from these elements:
        astro::SimpleOrbit orbit1(oe);    
   
        astro::TimeDelta dt(orbit1.getPeriod()/40.0);

        for(int i = 0; i <= 40; ++i)
        {
            // Get the state and print selected elements (position)
            auto state = orbit1.getState(et);
            std::cout << et.getETValue() << "\t" << state.r.x << "\t" << state.r.y << "\t" << state.r.z << std::endl;
            
            // Advance EphemerisTime with 1/40th period:
            et += dt;
 
        }
        
 
         
        
    }

    return 0;
}
