#include <iostream>
#include <limits>

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
            example = result["example"].as<int>();
        else
            std::cout << "Please give an example to run with -e <example>" << std::endl;
    }
    catch(cxxopts::OptionException)
    {
        std::cout << options.help({""}) << std::endl;
        return -2;
    }

    // **********************************************************
    // Example 1 - TIME
    // **********************************************************
    if(example == 1)
    {
        astro::Spice().loadKernel("../data/spice/lsk/naif0012.tls");

        astro::EphemerisTime et1;
        std::cout << "Date/Time in ISO: " << et1.toISOUTCString() << ", in Julian Day(UTC): " << et1.toJDUTCString() << std::endl;

        et1 = astro::EphemerisTime::fromString("2018 February 22, 20:04:00 UTC");
        std::cout << "Date/Time in ISO: " << et1.toISOUTCString() << ", in Julian Day(UTC): " << et1.toJDUTCString() << std::endl;

        et1 = astro::EphemerisTime::fromString("2017-12-24T17:00:00.12");
        std::cout << "Date/Time in ISO: " << et1.toISOUTCString() << ", in Julian Day(UTC): " << et1.toJDUTCString() << std::endl;

        et1 = astro::EphemerisTime::fromString("2451515.2981 JD");
        std::cout << "Date/Time in ISO: " << et1.toISOUTCString() << ", in Julian Day(UTC): " << et1.toJDUTCString() << std::endl;

        et1 = astro::EphemerisTime::fromJDUTC(2451545.0);
        std::cout << "Date/Time in ISO: " << et1.toISOUTCString() << ", in Julian Day(UTC): " << et1.toJDUTCString() << std::endl;
    }

    // **********************************************************
    // Example 2 - Simple Orbit & Orbit Elements
    // **********************************************************
    if(example == 2)
    {
        astro::PosState state;
        state.r = Vec3(-6045.0, -3490.0, 2500.0);  // [km]
        state.v = Vec3(-3.457, 6.618, 2.533);       // [km/s]

        astro::EphemerisTime et(0);
        double mu_earth = 398600.0;

        astro::OrbitElements oe = astro::OrbitElements::fromStateVector(state, et, mu_earth);
        astro::SimpleOrbit orbit1(oe);
        std::cout << orbit1.getOrbitElements() << std::endl;
    }

    // **********************************************************
    // Example 3 - Simple Orbit Plot
    // **********************************************************
    if(example == 3)
    {
        astro::PosState state;
        state.r = Vec3(-6045.0, -3490.0, 2500.0);  // [km]
        state.v = Vec3(-3.457, 6.618, 2.533);       // [km/s]

        astro::EphemerisTime et(0);
        double mu_earth = 398600.0;

        astro::OrbitElements oe = astro::OrbitElements::fromStateVector(state, et, mu_earth);
        astro::SimpleOrbit orbit1(oe);
        astro::TimeDelta dt(orbit1.getPeriod() / 40.0);

        for(int i = 0; i <= 40; ++i)
        {
            auto s = orbit1.getState(et);
            std::cout << et.getETValue() << "\t" << s.r.x << "\t" << s.r.y << "\t" << s.r.z << std::endl;
            et += dt;
        }
    }

    // **********************************************************
    // Example 4 - Orbit propagation
    // **********************************************************
    if(example == 4)
    {
        astro::PosState state;
        state.r = Vec3(7283.46, 0.0, 0.0);              // [km]
        state.v = Vec3(0.0, 58311.7 / 7283.46, 0.0);   // [km/s]

        astro::EphemerisTime et(0);
        double mu_earth = 398600.0;

        astro::OrbitElements oe = astro::OrbitElements::fromStateVector(state, et, mu_earth);
        astro::SimpleOrbit orbit1(oe);

        astro::ODE ode;
        astro::Attractor a = { Vec3(0.0), mu_earth };
        ode.addAttractor(a);

        astro::Propagator<astro::ODE, astro::RKF45> pr(ode);
        auto resv = pr.doSteps(state, et, et + astro::TimeDelta(orbit1.getPeriod()), astro::TimeDelta(1.0));

        for(auto res : resv)
        {
            auto oe2 = astro::OrbitElements::fromStateVector(res.s, res.et, mu_earth);
            std::cout << res.et.getETValue() << "\t" << res.dt_next.value
                      << "\t" << res.s.r.x << "\t" << res.s.r.y << "\t" << res.s.r.z
                      << "\t" << oe2.h << "\t" << oe2.a << std::endl;
        }
    }

    // **********************************************************
    // Example 5 - Hermite Spline interpolation
    // **********************************************************
    if(example == 5)
    {
        astro::PosState s1(Vec3(0,0,0),      Vec3(0,4,0));
        astro::EphemerisTime et1(1000.0);
        astro::PosState s2(Vec3(1000,2000,0), Vec3(2,0,0));
        astro::EphemerisTime et2(2000.0);
        astro::PosState sn;
        for(double t = 0.0; t < 1.0; t += 0.1)
        {
            astro::hermite(s1, et1, s2, et2, et1 + astro::TimeDelta((et2 - et1).value * t), sn);
            std::cout << sn.r.x << "\t" << sn.r.y << "\t" << sn.r.z << "\t"
                      << sn.v.x << "\t" << sn.v.y << "\t" << sn.v.z << std::endl;
        }
    }

    // **********************************************************
    // Example 6 - Integration of rotations
    // **********************************************************
    if(example == 6)
    {
        double w = 0.1;
        astro::RotState rs;
        rs.q = Quat(1.0, 0.0, 0.0, 0.0);  // identity (GLM: w first)
        rs.w = Vec3(w, 0.0, 0.0);           // rotation about global X at 0.1 rad/s

        astro::RotODE rode;
        astro::EphemerisTime et(12345);
        double T = astro::PIHALF / w;
        astro::EphemerisTime et2 = et + TimeDelta(T);
        astro::TimeDelta dt(1.0 / 60.0);

        auto resv = astro::PCDM::doSteps(rode, rs, et, et2, dt);

        // After a 90° rotation about X: x stays, y→z, z→-y
        Vec3 ex(1, 0, 0);
        Vec3 ey(0, 1, 0);
        Vec3 ez(0, 0, 1);

        const Quat& q = resv.back().rs.q;
        ex = q * ex;
        ey = q * ey;
        ez = q * ez;

        std::cout << ex << std::endl;
        std::cout << ey << std::endl;
        std::cout << ez << std::endl;
        std::cout << glm::length(q) << std::endl;
    }

    return 0;
}
