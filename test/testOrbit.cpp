// Reference from some of the tests:
// [1]  Orbital Mechanics for Engineering Students, 2nd Edition, Howard D. Curtis


#include "../astro/Orbit.cpp"
#include "../astro/State.h"
#include <gtest/gtest.h>

#include <cmath>

// Prettyprinters for vectors etc
#include "OrkExt.h"

using ork::vec3d;
using ork::mat3d;
using ork::mat4d;

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
};



OrbitTest::OrbitTest()
{

}

OrbitTest::~OrbitTest()
{

}

void OrbitTest::SetUp()
{
}

void OrbitTest::TearDown()
{
}

TEST_F(OrbitTest, OrbitalElementsFromState)
{
    double mu_earth = 398600.0;

    // [1], Example 4.3:
    astro::State   state;
    state.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
    state.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]

    astro::OrbitElements oe = astro::OrbitElements::fromStateVectorOE(state, mu_earth);
    // Compare with [1] (Those are quite rounded numbers, so we cant be too accurate)
    ASSERT_LT(fabs(oe.h- 58310), 2.0);
    ASSERT_LT(fabs(oe.i * astro::DEGPERRAD - 153.2), 1.0E-1);
    ASSERT_LT(fabs(oe.omega * astro::DEGPERRAD - 255.3), 1.0E-1);
    ASSERT_LT(fabs(oe.e - 0.1712), 1.0E-4);
    ASSERT_LT(fabs(oe.w * astro::DEGPERRAD - 20.07), 1.0E-2);
    ASSERT_LT(fabs(oe.theta * astro::DEGPERRAD - 28.45), 1.0E-2);

    // Invert the velocti, and see of we get a prograde orbit.
    double incl_old = oe.i * astro::DEGPERRAD;
    state.v *= -1.0;
    oe = astro::OrbitElements::fromStateVectorOE(state, mu_earth);
    double incl_new = oe.i * astro::DEGPERRAD;
    // Retrograde
    ASSERT_GT(incl_old, 90.0);
    ASSERT_LT(incl_old, 180.0);
    // Prograde
    ASSERT_GT(incl_new, 0.0);
    ASSERT_LT(incl_new, 90.0);
}

TEST_F(OrbitTest, OrbitalElementsFromState2)
{
    double mu_earth = 398600.0;

    // [1], Example 4.3:
    astro::State   state;
    state.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
    state.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]

    std::cout << "Input state:" << std::endl;
    print(state);

    astro::OrbitElements oe = astro::OrbitElements::fromStateVectorOE(state, mu_earth);
    astro::OrbitElements oe2 = astro::OrbitElements::fromStateVectorSpice(state, mu_earth);
    
    print(oe);   
    print(oe2);
}

TEST_F(OrbitTest, OrbitalElementsFromStateCornerCases)
{
    double mu = 398600.0;

    // [1], Example 4.3:
    astro::State   state;
    state.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
    state.v = vec3d(-3.457, 6.618, 2.533);    
    double v = state.v.length();
    double r = state.r.length();
    // Negative mass
    astro::OrbitElements oe;
	ASSERT_THROW(oe = astro::OrbitElements::fromStateVectorOE(state, -1000.0), astro::AstroException);
    ASSERT_THROW(oe = astro::OrbitElements::fromStateVectorSpice(state, -1000.0), astro::SpiceException);

    // Zero velocity
    state.v = vec3d::ZERO;
    ASSERT_THROW(oe = astro::OrbitElements::fromStateVectorOE(state, mu), astro::AstroException);
    ASSERT_THROW(oe = astro::OrbitElements::fromStateVectorSpice(state, mu), astro::SpiceException);

    // Zero radius
    state.v = vec3d(-3.457, 6.618, 2.533);
    state.r = vec3d::ZERO;
    ASSERT_THROW(oe = astro::OrbitElements::fromStateVectorOE(state, mu), astro::AstroException);
    ASSERT_THROW(oe = astro::OrbitElements::fromStateVectorSpice(state, mu), astro::SpiceException);

    // circular orbit:
    r = (6378.0 + 400);
    state.r = vec3d(r, 0.0, 0.0);
    v = 7.66895;
    state.v = vec3d(0, v/sqrt(2), v/sqrt(2));
    ASSERT_NO_THROW(oe = astro::OrbitElements::fromStateVectorOE(state, mu));
    print(oe);
    ASSERT_NO_THROW(oe = astro::OrbitElements::fromStateVectorSpice(state, mu));
    print(oe);

    
}

#if 0
TEST_F(OrbitTest, OrbitalElementsFromStateOpt)
{
    double mu_earth = 398600.0;

    // [1], Example 4.3:
    astro::State   state;
    state.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
    state.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]

    astro::OrbitElements oe = astro::OrbitElements::fromstateVectorOEOpt(state, mu_earth);
    // Compare with [1] (Those are quite rounded numbers, so we cant be too accurate)
    ASSERT_LT(fabs(oe.h- 58310), 2.0);
    ASSERT_LT(fabs(oe.i * astro::DEGPERRAD - 153.2), 1.0E-1);
    ASSERT_LT(fabs(oe.omega * astro::DEGPERRAD - 255.3), 1.0E-1);
    ASSERT_LT(fabs(oe.e - 0.1712), 1.0E-4);
    ASSERT_LT(fabs(oe.w * astro::DEGPERRAD - 20.07), 1.0E-2);
    ASSERT_LT(fabs(oe.theta * astro::DEGPERRAD - 28.45), 1.0E-2);

    // Invert the velocti, and see of we get a prograde orbit.
    double incl_old = oe.i * astro::DEGPERRAD;
    state.v *= -1.0;
    oe = astro::OrbitElements::fromstateVectorOEOpt(state, mu_earth);
    double incl_new = oe.i * astro::DEGPERRAD;
    // Retrograde
    ASSERT_GT(incl_old, 90.0);
    ASSERT_LT(incl_old, 180.0);
    // Prograde
    ASSERT_GT(incl_new, 0.0);
    ASSERT_LT(incl_new, 90.0);

}
#endif
// The benchmark timing for the first implementation
// Do 1 000 000 evaluations of the algorithm
TEST_F(OrbitTest, OrbitalElementsOEBenchmarkTime)
{
    double mu_earth = 398600.0;

    // [1], Example 4.3:
    astro::State   state;
    state.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
    state.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]

    astro::OrbitElements oe;
    for(auto i = 0; i < 1000000; ++i)
    {
        oe = astro::OrbitElements::fromStateVectorOE(state, mu_earth);
    }

}

TEST_F(OrbitTest, OrbitalElementsSpiceBenchmarkTime)
{
    double mu_earth = 398600.0;

    // [1], Example 4.3:
    astro::State   state;
    state.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
    state.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]

    astro::OrbitElements oe;
    for(auto i = 0; i < 1000000; ++i)
    {
        oe = astro::OrbitElements::fromStateVectorSpice(state, mu_earth);
    }

}

#if 0
TEST_F(OrbitTest, OrbitalElementsOptimizedTime)
{
    double mu_earth = 398600.0;

    // [1], Example 4.3:
    astro::State   state;
    state.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
    state.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]

    astro::OrbitElements oe;
    for(auto i = 0; i < 10000000; ++i)
    {
        oe = astro::OrbitElements::fromstateVectorOEOpt(state, mu_earth);
    }

}
#endif

TEST_F(OrbitTest, OrkRotationTest)
{
    vec3d I = vec3d::UNIT_X;
    vec3d J = vec3d::UNIT_Y;
    vec3d K = vec3d::UNIT_Z;

    // Rotate I 90 degrees about z axis should produce J
    // Note, ork::mat4 uses degrees!
    mat3d Q = mat4d::rotatez(90.0).mat3x3();
    vec3d rot = Q * I;
    ASSERT_LT((rot - J).length(), 1.0E-10);

    // Rotate Z -90 degrees about X should give J
    Q = mat4d::rotatex(-90.0).mat3x3();
    rot = Q * K;
    ASSERT_LT((rot - J).length(), 1.0E-10);

    // Rotate J 90 degrees about I shuold give K
    Q = mat4d::rotatex(90.0).mat3x3();
    rot = Q * J;
    ASSERT_LT((rot - K).length(), 1.0E-10);

}

TEST_F(OrbitTest, OrbitElementsToStateVector)
{
    double rpd = astro::RADPERDEG;
    double mu_earth = 398600.0;

    // [1], example 4.7:
    astro::OrbitElements oe;
    oe.h = 80000; // km²/s
    oe.e = 1.4;
    oe.i = 30.0 * rpd;
    oe.omega = 40.0 * rpd;
    oe.w = 60.0 * rpd;
    oe.theta = 30 * rpd;
    
    astro::State state = oe.toStateVectorOE(mu_earth);


    astro::OrbitElements oe2 = astro::OrbitElements::fromStateVectorOE(state, mu_earth);


    ASSERT_LT(fabs(oe2.h-oe.h), 1.0E-5);
    ASSERT_LT(fabs(oe2.i-oe.i), 1.0E-5);
    ASSERT_LT(fabs(oe2.omega-oe.omega), 1.0E-5);
    ASSERT_LT(fabs(oe2.w-oe.w), 1.0E-5);
    ASSERT_LT(fabs(oe2.theta-oe.theta), 1.0E-5);
    ASSERT_LT(fabs(oe2.e-oe.e), 1.0E-5);


}

TEST_F(OrbitTest, OrbitElementsToStateVector2)
{
    double mu_earth = 398600.0;

    // [1], Example 4.3:
    astro::State   state;
    state.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
    state.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]

    astro::OrbitElements oe = astro::OrbitElements::fromStateVectorOE(state, mu_earth);

    astro::State state2 = oe.toStateVectorOE(mu_earth);







}
