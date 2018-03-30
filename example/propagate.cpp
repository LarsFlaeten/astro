#include <iostream>
#include <string>
#include <iomanip>

#include <cxxopts.hpp>


#include <astro/SpiceCore.h>
#include <astro/Time.h>
#include <astro/State.h>
#include <astro/Orbit.h>
#include <astro/ODE.h>
#include <astro/Propagator.h>


int main(int argc, char **argv)
{
    cxxopts::Options options(argv[0], "Astrodynamics Library\n(c) 2018 Lars Flæten");
    options
      	.positional_help("RK1/RK2/RK3/RK4/RKF45/RKF78")
		.show_positional_help();

    std::string method;
    int periods;
    double DT;
    double tolerance;
    bool p_sta = false;
    bool p_oe = false;
    bool summary = false;
    try
    {
        options.add_options()
            ("h,help", "Print help")
            ("positional",
                "Positional arguments: these are the arguments that are entered "
                "without an option", cxxopts::value<std::vector<std::string>>())
            ("m,method", "Choose numerical method", cxxopts::value<std::string>())
            ("t,num_periods", "Number of Periods to integrate", cxxopts::value<int>(periods)->default_value("1"))
            ("d,delta_t", "Initial (for adaptive) or constant time step", cxxopts::value<double>(DT))
            ("s,summary", "Prints a summery of the method and parameters", cxxopts::value<bool>(summary)->default_value("false")->implicit_value("true"))
            ("p,print_state", "Prints the states of the result")
            ("l,tolerance", "Set tolerance for adaptive methods", cxxopts::value<double>(tolerance)->default_value("1.0E-8"))
            ("o,print_oel", "Prints the orbital elements of each result", cxxopts::value<bool>(p_oe)->default_value("false")->implicit_value("true"));
            
        options.parse_positional({"method", "positional"});
        auto result = options.parse(argc, argv);
        if(result.count("help"))
        {
            std::cout << options.help({""}) << std::endl;
            return 0;
        }

        if(result.count("method"))
        {
            method = result["method"].as<std::string>();
            
        } else
        {
            std::cout << "Please provide method (RKF45,....)" << std::endl;
            return -2;
        }
        if(result.count("print_state"))
            p_sta = true;
        
        if(result.count("delta_t"))
        {
            DT = result["delta_t"].as<double>();
            if(DT <= 0.0)
            {
                std::cout << "Negative or zero dt is not allowed" << std::endl;
                return -2;
                
            }
        }
        else
        {
            std::cout << "please supply initial/constant dt (-d / --delta_t)" << std::endl;
            return -2;
        }

            

        if(periods < 1)
        {
            std::cout << "Number of periods needs to be >1" << std::endl;
            return -2;
        }

    }
    catch(cxxopts::OptionException){
        std::cout << options.help({""}) << std::endl;
        return -2;
        
    }


    // **************************************************************************
    // Orbit propagation
    // This example shows how to propagate an orbit with numerical integration
    // **************************************************************************
   	// Roughly same orbit as 2&3, but with zero inclination,
	// and mean anomaly at epoch at periapsis
	astro::PosState   state0;
    state0.r = vec3d(7283.46, 0.0, 0.0);  //[km]
    state0.v = vec3d(0.0, 58311.7/7283.46, 0.0);      //[km/s]

    astro::EphemerisTime et0(0); // et = 0 represents the J2000 epoch
    double mu_earth = 398600.0;
		
	// Create a reference orbit for comparison
    astro::OrbitElements oe0 = astro::OrbitElements::fromStateVector(state0, et0, mu_earth);
    astro::SimpleOrbit orbit1(oe0);    
   	
    // The differential equation
    astro::ODE ode(mu_earth);

    // The time period to integrate	
    astro::TimeDelta period(orbit1.getPeriod() * periods);
    astro::EphemerisTime et1 = et0 + period;
 
    // A vector of orbital elements
    std::vector<astro::OrbitElements> oes;    
	
    if(method.compare("RKF45")==0)
    {
        astro::Propagator<astro::RKF45, astro::RKF45::Result> pr(ode);
        //astro::RKF45::setTolerance(tolerance);       

        auto resv = pr.doSteps(state0, et0, et1, astro::TimeDelta(DT));

        for( auto res : resv)
            oes.push_back(astro::OrbitElements::fromStateVector(res.s, res.et, mu_earth));    

        if(p_sta)
            for(astro::RKF45::Result res : resv)
                std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << res.et.getETValue() << "\t" << res.dt_next.value << "\t" <<
                res.s.r.x << "\t" << res.s.r.y << "\t" << res.s.r.z << "\t" <<
                res.s.v.x << "\t" << res.s.v.y << "\t" << res.s.v.z << std::endl;
    }
    else if(method.compare("RKF78")==0)
    {

        astro::Propagator<astro::RKF78, astro::RKF78::Result> pr(ode);
        //astro::RKF78::setTolerance(tolerance);

        auto resv = pr.doSteps(state0, et0, et1, astro::TimeDelta(DT));

        for( auto res : resv)
            oes.push_back(astro::OrbitElements::fromStateVector(res.s, res.et, mu_earth));    

        if(p_sta)
            for(astro::RKF78::Result res : resv)
                std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << res.et.getETValue() << "\t" << res.dt_next.value << "\t" <<
                res.s.r.x << "\t" << res.s.r.y << "\t" << res.s.r.z << "\t" <<
                res.s.v.x << "\t" << res.s.v.y << "\t" << res.s.v.z << std::endl;

   }     
    else if(method.compare("RK1")==0)
    {
        astro::Propagator<astro::RK<1>, astro::RK<1>::Result> pr(ode);
       
        auto resv = pr.doSteps(state0, et0, et1, astro::TimeDelta(DT));

        for( auto res : resv)
            oes.push_back(astro::OrbitElements::fromStateVector(res.s, res.et, mu_earth));    

        if(p_sta)
            for(auto res : resv)
                std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << res.et.getETValue() << "\t" << 
                res.s.r.x << "\t" << res.s.r.y << "\t" << res.s.r.z << "\t" <<
                res.s.v.x << "\t" << res.s.v.y << "\t" << res.s.v.z << std::endl;
    } 
    else if(method.compare("RK2")==0)
    {
        astro::Propagator<astro::RK<2>, astro::RK<2>::Result> pr(ode);
       
        auto resv = pr.doSteps(state0, et0, et1, astro::TimeDelta(DT));

        for( auto res : resv)
            oes.push_back(astro::OrbitElements::fromStateVector(res.s, res.et, mu_earth));    

        if(p_sta)
            for(auto res : resv)
                std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << res.et.getETValue() << "\t" << 
                res.s.r.x << "\t" << res.s.r.y << "\t" << res.s.r.z << "\t" <<
                res.s.v.x << "\t" << res.s.v.y << "\t" << res.s.v.z << std::endl;
    } 
    else if(method.compare("RK3")==0)
    {
        astro::Propagator<astro::RK<3>, astro::RK<3>::Result> pr(ode);
       
        auto resv = pr.doSteps(state0, et0, et1, astro::TimeDelta(DT));

        for( auto res : resv)
            oes.push_back(astro::OrbitElements::fromStateVector(res.s, res.et, mu_earth));    

        if(p_sta)
            for(auto res : resv)
                std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << res.et.getETValue() << "\t" << 
                res.s.r.x << "\t" << res.s.r.y << "\t" << res.s.r.z << "\t" <<
                res.s.v.x << "\t" << res.s.v.y << "\t" << res.s.v.z << std::endl;
    } 
    else if(method.compare("RK4")==0)
    {
        astro::Propagator<astro::RK<4>, astro::RK<4>::Result> pr(ode);
       
        auto resv = pr.doSteps(state0, et0, et1, astro::TimeDelta(DT));

        for( auto res : resv)
            oes.push_back(astro::OrbitElements::fromStateVector(res.s, res.et, mu_earth));    

        if(p_sta)
            for(auto res : resv)
                std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << res.et.getETValue() << "\t" << 
                res.s.r.x << "\t" << res.s.r.y << "\t" << res.s.r.z << "\t" <<
                res.s.v.x << "\t" << res.s.v.y << "\t" << res.s.v.z << std::endl;
    } 
    // Test of ODEints implementation.
    else if(method.compare("RK4O")==0)
    {
        runge_kutta4< astro::PosState, double, astro::PosState, double, vector_space_algebra > rk4;
        astro::PosState s = state0;
        if(p_sta)
            std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << et0.getETValue() << "\t" << 
            s.r.x << "\t" << s.r.y << "\t" << s.r.z << "\t" <<
            s.v.x << "\t" << s.v.y << "\t" << s.v.z << std::endl;
        
        astro::EphemerisTime et = et0;
        while(et < et1)
        {
            DT = std::min(DT, (et1 - et).value);
            
            rk4.do_step(ode, s, et.getETValue(), DT);
            et += DT;
    
            if(p_sta)
                std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << et.getETValue() << "\t" << 
                s.r.x << "\t" << s.r.y << "\t" << s.r.z << "\t" <<
                s.v.x << "\t" << s.v.y << "\t" << s.v.z << std::endl;
        } 
        //auto resv = pr.doSteps(state0, et0, et1, astro::TimeDelta(DT));

        //for( auto res : resv)
        //    oes.push_back(astro::OrbitElements::fromStateVector(res.s, res.et, mu_earth));    

        //if(p_sta)
        //    for(auto res : resv)
        //        std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << res.et.getETValue() << "\t" << 
        //        res.s.r.x << "\t" << res.s.r.y << "\t" << res.s.r.z << "\t" <<
        //        res.s.v.x << "\t" << res.s.v.y << "\t" << res.s.v.z << std::endl;
    } 
    else
    {
        std::cout << "Unkown method. Valid methods are RKF45" << std::endl;
        return -2;
    }
    
    if(p_oe)
    {
       std::cout << "oe.h" << "\t" << "oe.a" << "\t" << "oe.e" << "\t" <<
            "oe.i" << "\t" << "oe.omega" << "\t" << "oe.w" << "\t" << "oe.M0" << std::endl;

       for(auto oe : oes)
            std::cout << oe.h << "\t" << oe.a << "\t" << oe.e << "\t" <<
            oe.i << "\t" << oe.omega << "\t" << oe.w << "\t" << oe.M0 << std::endl;

    }


    if(summary)
    {
        std::cout << "Using method=\"" << method << "\"" << std::endl;
        std::cout << "Number of periods=" << periods << std::endl;
    }
    
    return 0;
}
