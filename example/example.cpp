#include <iostream>

#include <cxxopts.hpp>

#include <astro/Util.h>
#include <astro/SpiceCore.h>
#include <astro/Time.h>
#include <astro/State.h>
#include <astro/Orbit.h>
#include <astro/ODE.h>
#include <astro/Propagator.h>
#include <astro/PCDM.h>
#include <astro/Interpolate.h>

using namespace astro;

int main(int argc, char **argv)
{
    int example = 0;
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
        astro::PosState   state;
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
        astro::PosState   state;
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

    // **************************************************************************
    // Example 4 - Orbit propagation
    // This example shows how to propagate an orbit with numerical integration
    // **************************************************************************
    if(example == 4)
    {
		// Roughly same orbit as 2&3, but with zero inclination,
		// and mean anomaly at epoch at periapsis
		astro::PosState   state;
        state.r = vec3d(7283.46, 0.0, 0.0);  //[km]
        state.v = vec3d(0.0, 58311.7/7283.46, 0.0);      //[km/s]

        astro::EphemerisTime et(0); // et = 0 represents the J2000 epoch
        double mu_earth = 398600.0;
		
		// Create a reference orbit for comparison
        astro::OrbitElements oe = astro::OrbitElements::fromStateVector(state, et, mu_earth);
        astro::SimpleOrbit orbit1(oe);    
        //std::cout << orbit1.getOrbitElements() << std::endl;
   		
		astro::ODE ode;
	    astro::Attractor a = {vec3d::ZERO, mu_earth};
        ode.addAttractor(a);

        astro::Propagator<astro::ODE, astro::RKF45> pr(ode);
        astro::TimeDelta period(orbit1.getPeriod());
        astro::TimeDelta dt(1.0);

        // Do the integration:      
        auto resv = pr.doSteps(state, et, et+period, dt);

        // Dump positions to stdout:
        for(auto res : resv)
        {
            std::cout << res.et.getETValue() << "\t" << res.dt_next.value << "\t" << res.s.r.x << "\t" << res.s.r.y << "\t" << res.s.r.z ;

            // Angular momentum and semimajor axis shoud be stable:
            auto oe = astro::OrbitElements::fromStateVector(res.s, res.et, mu_earth);

            std::cout << "\t" << oe.h << "\t" << oe.a << std::endl;
        }


	}

    // **************************************************************************
    // Example 5 - Hermite Spline interpolation
    // This example shows how to interpolate between known states
    // **************************************************************************
    if(example == 5)
    {
        astro::PosState s1 = {vec3d(0,0,0), vec3d(0,4,0)};
        astro::EphemerisTime et1(1000.0);
        astro::PosState s2 = {vec3d(1000,2000,0), vec3d(2,0,0)};	        
        astro::EphemerisTime et2(2000.0);
        astro::PosState sn;
        for(double t = 0.0; t < 1.0; t += 0.1)
        {
            astro::hermite(s1, et1, s2, et2, et1 + astro::TimeDelta((et2-et1).value*t), sn);
            std::cout << sn.r.x << "\t" << sn.r.y << "\t" << sn.r.z << "\t";
            std::cout << sn.v.x << "\t" << sn.v.y << "\t" << sn.v.z << std::endl;
        }
    }

    // **************************************************************************
    // Example 6 - Integration rotations
    // This example shows how to numerically integrate rotation states
    // **************************************************************************
    if(example == 6)
    {
        double w = 0.1;
        astro::RotState rs;
        rs.q = quatd(0, 0, 0, 1); // "Unit" quaternion
        rs.w = vec3d(w, 0.0, 0.0); // rotation about global X by 0.1 rad/s

        // The differential equation for rotations:
        astro::RotODE rode;

        astro::EphemerisTime et(12345);
        // Find the time when we shold be rotated 90 degrees by x:
        double T = astro::PIHALF / w;
        astro::EphemerisTime et2 = et + TimeDelta(T);
        

        astro::TimeDelta dt(1/60.0);
        //astro::TimeDelta dt(1.0/60);

        
        auto resv = astro::PCDM::doSteps(rode, rs, et, et2, dt);

        // After this, x will be x, y will be z, and z will be -y.
        vec3d ex = vec3d::UNIT_X;
        vec3d ey = vec3d::UNIT_Y;
        vec3d ez = vec3d::UNIT_Z;

        ex = resv.back().rs.q * ex;
        ey = resv.back().rs.q * ey;
        ez = resv.back().rs.q * ez;

        std::cout << ex << std::endl;
        std::cout << ey << std::endl;
        std::cout << ez << std::endl;

        std::cout << resv.back().rs.q.length() << std::endl;



       
    }


    
    return 0;
}
