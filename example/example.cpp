#include <iostream>

#include <cxxopts.hpp>


#include <astro/SpiceCore.h>
#include <astro/Time.h>
#include <astro/State.h>
#include <astro/Orbit.h>
#include <astro/ODE.h>
#include <astro/Propagator.h>
#include <astro/Interpolate.h>

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
   		
		astro::ODE ode(mu_earth);
		
        astro::Propagator<astro::RKF45, astro::RKF45::Result> pr(ode);
        astro::TimeDelta period(orbit1.getPeriod());
        astro::TimeDelta dt(1.0);
        
       
        // create vector of results, and start with initial condition: 
        std::vector<astro::RKF45::Result> results;
        results.push_back({state, et, dt}); 
       
        // creat vector of orbit properties, so we can track those over time: 
        std::vector<astro::OrbitElements>  oes;
        oes.push_back(oe);
	    

 	    astro::EphemerisTime eti = et;
        astro::PosState statei = state;
        int i = 0;
        while(1)
        {
            ++i;
            auto res = pr.doStep(statei, eti, dt);
		   
            oes.push_back(astro::OrbitElements::fromStateVector(statei, eti, mu_earth));
            results.push_back(res);
             
            dt = res.dt_next;
            eti = res.et;
            statei = res.s;
            
             
            if(res.et > et + period)
                break;
        }

        // Dump positions to stdout:
#if 0
        for(astro::RKF45::Result res : results)
        {
            std::cout << res.et.getETValue() << "\t" << res.s.r.x << "\t" << res.s.r.y << "\t" << res.s.r.z << std::endl;
        }
#endif

        // Dump angular momentum and semimajor axis to stdout:
        for(astro::OrbitElements oen : oes)
            std::cout << oen.h << "\t" << oen.a << std::endl;
	}

    // **************************************************************************
    // Example 5 - Hermite Spline interpolation
    // This example shows how to interpolate between known states
    // **************************************************************************
    if(example == 5)
    {
        astro::PosState s1 = {vec3d(0,0,0), vec3d(0,2,0)};
        astro::PosState s2 = {vec3d(1,2,0), vec3d(1,0,0)};	        
        astro::PosState sn;
        for(double t = 0.0; t < 1.0; t += 0.1)
        {
            astro::hermite(s1, s2, t, sn);
            std::cout << sn.r.x << "\t" << sn.r.y << "\t" << sn.r.z << "\t";
            std::cout << sn.v.x << "\t" << sn.v.y << "\t" << sn.v.z << std::endl;





        }
    }
    return 0;
}
