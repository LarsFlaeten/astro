#include "../astro/State.h"
#include "../astro/Time.h"
#include "../astro/Propagator.h"
#include "../astro/RKF78.cpp"
#include "../astro/RKF45.cpp"
#include "../astro/Orbit.h"
#include "../astro/ODE.h"

#include <gtest/gtest.h>

#include <cmath>

using namespace boost::numeric::odeint;

class NumIntTest : public ::testing::Test {

protected:
    NumIntTest();

    virtual ~NumIntTest();

    // Code here will be called immediately after the constructor (right
    // before each test).
    virtual void SetUp();

    // Code here will be called immediately after each test (right
    // before the destructor).
    virtual void TearDown();

    astro::State    state0;
    astro::OrbitElements oe0;
    astro::EphemerisTime et0;
    double mu_earth;
    astro::ODE  ode0;
};



NumIntTest::NumIntTest()
  :  et0(0), mu_earth(398600.0), ode0(mu_earth)
{
    
    state0.r = vec3d(7283.46, 0.0, 0.0);  //[km]
    state0.v = vec3d(0.0, 58311.7/7283.46, 0.0);      //[km/s]

    oe0 = astro::OrbitElements::fromStateVector(state0, et0, mu_earth);

}

NumIntTest::~NumIntTest()
{

}

void NumIntTest::SetUp()
{
}

void NumIntTest::TearDown()
{
}

TEST_F(NumIntTest, RKF78BenchMarkTest)
{

    astro::SimpleOrbit orbit1(oe0);    
    
    double DT = 0.1;

    // The time period to integrate	
    astro::TimeDelta period(orbit1.getPeriod() * 5000);
    astro::EphemerisTime et1 = et0 + period;
    
    astro::Propagator<astro::RKF78, astro::RKF78::Result> pr(ode0);
    astro::RKF78::setTolerance(1.0E-8);       

    auto resv = pr.doSteps(state0, et0, et1, astro::TimeDelta(DT));


    	

}

TEST_F(NumIntTest, RKF45BenchMarkTest)
{

    astro::SimpleOrbit orbit1(oe0);    
    
    double DT = 0.1;

    // The time period to integrate	
    astro::TimeDelta period(orbit1.getPeriod() * 5000);
    astro::EphemerisTime et1 = et0 + period;
    
    astro::Propagator<astro::RKF45, astro::RKF45::Result> pr(ode0);
    astro::RKF45::setTolerance(1.0E-8);       

    auto resv = pr.doSteps(state0, et0, et1, astro::TimeDelta(DT));


    	

}



