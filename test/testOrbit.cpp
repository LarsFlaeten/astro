// Reference from some of the tests:
// [1]  Orbital Mechanics for Engineering Students, 2nd Edition, Howard D. Curtis


#include "../astro/Orbit.cpp"
#include "../astro/State.h"
#include <gtest/gtest.h>

#include <cmath>

// Prettyprinters for vectors etc
//#include "OrkExt.h"

using mork::vec3d;
using mork::mat3d;
using mork::mat4d;

using namespace astro;

class OrbitTest : public ::testing::Test {

protected:
    OrbitTest();

    virtual ~OrbitTest();

    // Code here will be called immediately after the constructor (right
    // before each test).
    virtual void SetUp();

    // Code here will be called immediately after each test (right
    // before the destructor).
    virtual void TearDown();

    EphemerisTime et;
    double mu_earth;
    astro::OrbitElements    oe_ell;
    astro::OrbitElements    oe_hyp;

};



OrbitTest::OrbitTest()
 : et(), mu_earth(398600.0)
{
       
}

OrbitTest::~OrbitTest()
{

}

void OrbitTest::SetUp()
{
    // [1], Example 4.3:
    astro::PosState   state;
    state.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
    state.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]

    oe_ell = astro::OrbitElements::fromStateVectorOE(state, et, mu_earth);
    
    // [1], example 4.7:
    double rpd = astro::RADPERDEG;
    astro::OrbitElements oe;
    oe.h = 80000; // km²/s
    oe.e = 1.4;
    oe.i = 30.0 * rpd;
    oe.omega = 40.0 * rpd;
    oe.w = 60.0 * rpd;
    oe.M0 = astro::OrbitElements::meanAnomalyFromTrueAnomaly(30 * rpd, oe.e);
    oe.epoch = et;
    oe.mu = mu_earth;
    
    oe_hyp = oe;
    oe_hyp.computeDerivedQuantities();


}

void OrbitTest::TearDown()
{
}

TEST_F(OrbitTest, InitTest)
{
    astro::SimpleOrbit o(oe_ell);

    ASSERT_GT(o.getPeriod(), 0.0);
}

TEST_F(OrbitTest, GetPeriodTest)
{
    astro::SimpleOrbit o1(oe_ell);
    astro::SimpleOrbit o2(oe_hyp);

    ASSERT_LT(fabs(o1.getPeriod()-8198.86), 1.0E-2);

    ASSERT_LT(o2.getPeriod(), 0.0);
    
}

TEST_F(OrbitTest, IsPeriodicTest)
{
    astro::SimpleOrbit o1(oe_ell);
    astro::SimpleOrbit o2(oe_hyp);

    ASSERT_EQ(o1.isPeriodic(), true);
    ASSERT_EQ(o2.isPeriodic(), false);


    
}




