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


    EphemerisTime et;
    double mu_earth;
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
}

void OrbitTest::TearDown()
{
}

TEST_F(OrbitTest, OrbitalElementsFromState)
{

    // [1], Example 4.3:
    astro::State   state;
    state.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
    state.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]

    astro::OrbitElements oe = astro::OrbitElements::fromStateVectorOE(state, et, mu_earth);
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
    oe = astro::OrbitElements::fromStateVectorOE(state, et, mu_earth);
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

    // [1], Example 4.3:
    astro::State   state;
    state.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
    state.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]

    //std::cout << "Input state:" << std::endl;
    //print(state);

    astro::OrbitElements oe = astro::OrbitElements::fromStateVectorOE(state, et, mu_earth);
    astro::OrbitElements oe2 = astro::OrbitElements::fromStateVectorSpice(state, et, mu_earth);
    ASSERT_LT(fabs(oe.h-oe2.h), 0.01);
    ASSERT_LT(fabs(oe.i-oe2.i), 0.001);
    ASSERT_LT(fabs(oe.omega-oe2.omega), 0.001);
    ASSERT_LT(fabs(oe.e-oe2.e), 1.0E-6);
    ASSERT_LT(fabs(oe.w-oe2.w), 0.0001);
    ASSERT_LT(fabs(oe.theta-oe2.theta), 0.0001);
    ASSERT_LT(fabs(oe.M0-oe2.M0), 0.0001);
    ASSERT_LT(fabs(oe.rp-oe2.rp), 0.0001);
    
    //print(oe);   
    //print(oe2);
}

TEST_F(OrbitTest, OrbitalElementsFromStateCornerCases)
{
	double dpr = astro::DEGPERRAD;
    double mu = mu_earth;

    // [1], Example 4.3:
    astro::State   state;
    state.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
    state.v = vec3d(-3.457, 6.618, 2.533);    
    double v = state.v.length();
    double r = state.r.length();
    // Negative mass
    astro::OrbitElements oe;
	ASSERT_THROW(oe = astro::OrbitElements::fromStateVectorOE(state, et, -1000.0), astro::AstroException);
    ASSERT_THROW(oe = astro::OrbitElements::fromStateVectorSpice(state, et, -1000.0), astro::SpiceException);

    // Zero velocity
    state.v = vec3d::ZERO;
    ASSERT_THROW(oe = astro::OrbitElements::fromStateVectorOE(state, et, mu), astro::AstroException);
    ASSERT_THROW(oe = astro::OrbitElements::fromStateVectorSpice(state, et, mu), astro::SpiceException);

    // Zero radius
    state.v = vec3d(-3.457, 6.618, 2.533);
    state.r = vec3d::ZERO;
    ASSERT_THROW(oe = astro::OrbitElements::fromStateVectorOE(state, et, mu), astro::AstroException);
    ASSERT_THROW(oe = astro::OrbitElements::fromStateVectorSpice(state, et, mu), astro::SpiceException);

    // circular orbit, incl = 45:
    astro::OrbitElements oe2;
    r = (6378.0 + 400);
    state.r = vec3d(r, 0.0, 0.0);
    v = 7.66895;
    state.v = vec3d(0, v/sqrt(2), v/sqrt(2));
    ASSERT_NO_THROW(oe = astro::OrbitElements::fromStateVectorOE(state, et, mu));
    //print(oe);
    ASSERT_NO_THROW(oe2 = astro::OrbitElements::fromStateVectorSpice(state, et, mu));
    //print(oe2);
    ASSERT_LT(fabs(oe.h-oe2.h), 0.01);
    ASSERT_LT(fabs(oe.i-oe2.i), 0.001);
    ASSERT_LT(fabs(oe.omega-oe2.omega), 0.001);
    ASSERT_LT(fabs(oe.e-oe2.e), 1.0E-6);
    ASSERT_LT(fabs(oe.w-oe2.w), 0.0001);
    ASSERT_LT(fabs(oe.theta-oe2.theta), 0.0001);
    ASSERT_LT(fabs(oe.M0-oe2.M0), 0.0001);
    ASSERT_LT(fabs(oe.rp-oe2.rp), 0.0001);
    
	// Elliptic orbit, incl = 45:
    r = (6378.0 + 400);
    state.r = vec3d(r, 0.0, 0.0);
    v = 7.66895;
    state.v = vec3d(0, v*1.1/sqrt(2), v*1.1/sqrt(2));
    ASSERT_NO_THROW(oe = astro::OrbitElements::fromStateVectorOE(state, et, mu));
    //print(oe);
    ASSERT_NO_THROW(oe2 = astro::OrbitElements::fromStateVectorSpice(state, et, mu));
    //print(oe2);
    ASSERT_LT(fabs(oe.h-oe2.h), 0.01);
    ASSERT_LT(fabs(oe.i-oe2.i), 0.001);
    ASSERT_LT(fabs(oe.omega-oe2.omega), 0.001);
    ASSERT_LT(fabs(oe.e-oe2.e), 1.0E-6);
    ASSERT_LT(fabs(oe.w-oe2.w), 0.0001);
    ASSERT_LT(fabs(oe.theta-oe2.theta), 0.0001);
    ASSERT_LT(fabs(oe.M0-oe2.M0), 0.0001);
    ASSERT_LT(fabs(oe.rp-oe2.rp), 0.0001);

    // hyperbolic orbit, incl = 45:
    r = (6378.0 + 400);
    state.r = vec3d(r, 0.0, 0.0);
    v = 7.66895;
    state.v = vec3d(0, v*1.5/sqrt(2), v*1.5/sqrt(2));
    ASSERT_NO_THROW(oe = astro::OrbitElements::fromStateVectorOE(state, et, mu));
    //print(oe);
    ASSERT_NO_THROW(oe2 = astro::OrbitElements::fromStateVectorSpice(state, et, mu));
    //print(oe2);
    ASSERT_LT(fabs(oe.h-oe2.h), 0.01);
    ASSERT_LT(fabs(oe.i-oe2.i), 0.001);
    ASSERT_LT(fabs(oe.omega-oe2.omega), 0.001);
    ASSERT_LT(fabs(oe.e-oe2.e), 1.0E-6);
    ASSERT_LT(fabs(oe.w-oe2.w), 0.0001);
    ASSERT_LT(fabs(oe.theta-oe2.theta), 0.0001);
    ASSERT_LT(fabs(oe.M0-oe2.M0), 0.0001);
    ASSERT_LT(fabs(oe.rp-oe2.rp), 0.0001);

    // hyperbolic orbit 2, anomaly != 0:
    // From [1], Ex. 4.7:
    state.r = vec3d(-4040, 4815, 3629);
    state.v = vec3d(-10.39, -4.772, 1.744);
    ASSERT_NO_THROW(oe = astro::OrbitElements::fromStateVectorOE(state, et, mu));
    //print(oe);
    ASSERT_NO_THROW(oe2 = astro::OrbitElements::fromStateVectorSpice(state, et, mu));
    //print(oe2);
	ASSERT_LT(fabs(oe.h-oe2.h), 0.01);
    ASSERT_LT(fabs(oe.i-oe2.i), 0.001);
    ASSERT_LT(fabs(oe.omega-oe2.omega), 0.001);
    ASSERT_LT(fabs(oe.e-oe2.e), 1.0E-6);
    ASSERT_LT(fabs(oe.w-oe2.w), 0.0001);
    ASSERT_LT(fabs(oe.theta-oe2.theta), 0.0001);
    ASSERT_LT(fabs(oe.M0-oe2.M0), 0.0001);
    ASSERT_LT(fabs(oe.rp-oe2.rp), 0.0001);
	ASSERT_LT(fabs(oe.h-80000), 30);
    ASSERT_LT(fabs(oe.i*dpr-30), 0.01);
    ASSERT_LT(fabs(oe.omega*dpr-40), 0.022);
    ASSERT_LT(fabs(oe.e-1.4), 1.0E-2);
    ASSERT_LT(fabs(oe.w*dpr-60), 0.02);
    ASSERT_LT(fabs(oe.theta*dpr-30), 0.01);






	// Equatorial circular orbit
	state.v = vec3d(0, v, 0);
	ASSERT_NO_THROW(oe = astro::OrbitElements::fromStateVectorOE(state, et, mu));
	ASSERT_NO_THROW(oe2 = astro::OrbitElements::fromStateVectorSpice(state, et, mu));
	//print(oe);
	//print(oe2);
    ASSERT_LT(fabs(oe.h-oe2.h), 0.01);
    ASSERT_LT(fabs(oe.i-oe2.i), 0.001);
    ASSERT_LT(fabs(oe.omega-oe2.omega), 0.001);
    ASSERT_LT(fabs(oe.e-oe2.e), 1.0E-6);
    ASSERT_LT(fabs(oe.w-oe2.w), 0.0001);
    ASSERT_LT(fabs(oe.theta-oe2.theta), 0.0001);
    ASSERT_LT(fabs(oe.M0-oe2.M0), 0.0001);
    ASSERT_LT(fabs(oe.rp-oe2.rp), 0.0001);

	// Equatorial elliptic orbit
	state.v = vec3d(0, v*1.1, 0);
    ASSERT_NO_THROW(oe = astro::OrbitElements::fromStateVectorOE(state, et, mu));
    ASSERT_NO_THROW(oe2 = astro::OrbitElements::fromStateVectorSpice(state, et, mu));
    //print(oe);
    //print(oe2);
    ASSERT_LT(fabs(oe.h-oe2.h), 0.01);
    ASSERT_LT(fabs(oe.i-oe2.i), 0.001);
    ASSERT_LT(fabs(oe.omega-oe2.omega), 0.001);
    ASSERT_LT(fabs(oe.e-oe2.e), 1.0E-6);
    ASSERT_LT(fabs(oe.w-oe2.w), 0.0001);
    ASSERT_LT(fabs(oe.theta-oe2.theta), 0.0001);
    ASSERT_LT(fabs(oe.M0-oe2.M0), 0.0001);
    ASSERT_LT(fabs(oe.rp-oe2.rp), 0.0001);

    // Equatorial hyperbolic orbit
    state.v = vec3d(0, v*1.5, 0);
    ASSERT_NO_THROW(oe = astro::OrbitElements::fromStateVectorOE(state, et, mu));
    ASSERT_NO_THROW(oe2 = astro::OrbitElements::fromStateVectorSpice(state, et, mu));
    //print(oe);
    //print(oe2);
    ASSERT_LT(fabs(oe.h-oe2.h), 0.01);
    ASSERT_LT(fabs(oe.i-oe2.i), 0.001);
    ASSERT_LT(fabs(oe.omega-oe2.omega), 0.001);
    ASSERT_LT(fabs(oe.e-oe2.e), 1.0E-6);
    ASSERT_LT(fabs(oe.w-oe2.w), 0.0001);
    ASSERT_LT(fabs(oe.theta-oe2.theta), 0.0001);
    ASSERT_LT(fabs(oe.M0-oe2.M0), 0.0001);
    ASSERT_LT(fabs(oe.rp-oe2.rp), 0.0001);


	// Polar circular orbit
	state.v = vec3d(0, 0, v);
    ASSERT_NO_THROW(oe = astro::OrbitElements::fromStateVectorOE(state, et, mu));
    ASSERT_NO_THROW(oe2 = astro::OrbitElements::fromStateVectorSpice(state, et, mu));
	//print(oe);
    //print(oe2);
	ASSERT_LT(fabs(oe.h-oe2.h), 0.01);
    ASSERT_LT(fabs(oe.i-oe2.i), 0.001);
    ASSERT_LT(fabs(oe.omega-oe2.omega), 0.001);
    ASSERT_LT(fabs(oe.e-oe2.e), 1.0E-6);
    ASSERT_LT(fabs(oe.w-oe2.w), 0.0001);
    ASSERT_LT(fabs(oe.theta-oe2.theta), 0.0001);
    ASSERT_LT(fabs(oe.M0-oe2.M0), 0.0001);
    ASSERT_LT(fabs(oe.rp-oe2.rp), 0.0001);

	// Polar elliptic orbit
	state.v = vec3d(0, 0, v*1.1);
    ASSERT_NO_THROW(oe = astro::OrbitElements::fromStateVectorOE(state, et, mu));
    ASSERT_NO_THROW(oe2 = astro::OrbitElements::fromStateVectorSpice(state, et, mu));
    //print(oe);
    //print(oe2);
    ASSERT_LT(fabs(oe.h-oe2.h), 0.01);
    ASSERT_LT(fabs(oe.i-oe2.i), 0.001);
    ASSERT_LT(fabs(oe.omega-oe2.omega), 0.001);
    ASSERT_LT(fabs(oe.e-oe2.e), 1.0E-6);
    ASSERT_LT(fabs(oe.w-oe2.w), 0.0001);
    ASSERT_LT(fabs(oe.theta-oe2.theta), 0.0001);
    ASSERT_LT(fabs(oe.M0-oe2.M0), 0.0001);
    ASSERT_LT(fabs(oe.rp-oe2.rp), 0.0001);

    // Polar hyperbolic orbit
    state.v = vec3d(0, 0, v*1.5);
    ASSERT_NO_THROW(oe = astro::OrbitElements::fromStateVectorOE(state, et, mu));
    ASSERT_NO_THROW(oe2 = astro::OrbitElements::fromStateVectorSpice(state, et, mu));
    //print(oe);
    //print(oe2);
    ASSERT_LT(fabs(oe.h-oe2.h), 0.01);
    ASSERT_LT(fabs(oe.i-oe2.i), 0.001);
    ASSERT_LT(fabs(oe.omega-oe2.omega), 0.001);
    ASSERT_LT(fabs(oe.e-oe2.e), 1.0E-6);
    ASSERT_LT(fabs(oe.w-oe2.w), 0.0001);
    ASSERT_LT(fabs(oe.theta-oe2.theta), 0.0001);
    ASSERT_LT(fabs(oe.M0-oe2.M0), 0.0001);
    ASSERT_LT(fabs(oe.rp-oe2.rp), 0.0001);





    
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
    // [1], Example 4.3:
    astro::State   state;
    state.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
    state.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]

    astro::OrbitElements oe;
    for(auto i = 0; i < 1000000; ++i)
    {
        oe = astro::OrbitElements::fromStateVectorOE(state, et, mu_earth);
    }

}

TEST_F(OrbitTest, OrbitalElementsSpiceBenchmarkTime)
{

    // [1], Example 4.3:
    astro::State   state;
    state.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
    state.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]

    astro::OrbitElements oe;
    for(auto i = 0; i < 1000000; ++i)
    {
        oe = astro::OrbitElements::fromStateVectorSpice(state, et, mu_earth);
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

    // [1], example 4.7:
    astro::OrbitElements oe;
    oe.h = 80000; // km²/s
    oe.e = 1.4;
    oe.i = 30.0 * rpd;
    oe.omega = 40.0 * rpd;
    oe.w = 60.0 * rpd;
    oe.theta = 30 * rpd;
    oe.epoch = et;
    oe.mu = mu_earth;
    oe.rp = oe.h*oe.h / oe.mu * (1.0 / (1.0 + oe.e));

    astro::State state = oe.toStateVectorOE();


    astro::OrbitElements oe2 = astro::OrbitElements::fromStateVectorOE(state, et, mu_earth);


    ASSERT_LT(fabs(oe2.h-oe.h), 1.0E-5);
    ASSERT_LT(fabs(oe2.i-oe.i), 1.0E-5);
    ASSERT_LT(fabs(oe2.omega-oe.omega), 1.0E-5);
    ASSERT_LT(fabs(oe2.w-oe.w), 1.0E-5);
    ASSERT_LT(fabs(oe2.theta-oe.theta), 1.0E-5);
    ASSERT_LT(fabs(oe2.e-oe.e), 1.0E-5);


}

TEST_F(OrbitTest, OrbitElementsToStateVector2)
{
    // [1], Example 4.3:
    astro::State   state;
    state.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
    state.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]

    astro::OrbitElements oe = astro::OrbitElements::fromStateVectorOE(state, et, mu_earth);

    astro::State state2 = oe.toStateVectorOE();
    //print(state);
    //print(state2);
	ASSERT_LT(fabs(state.r.x-state2.r.x), 1.0E-5);
    ASSERT_LT(fabs(state.r.y-state2.r.y), 1.0E-5);
    ASSERT_LT(fabs(state.r.z-state2.r.z), 1.0E-5);
    ASSERT_LT(fabs(state.v.x-state2.v.x), 1.0E-5);
    ASSERT_LT(fabs(state.v.y-state2.v.y), 1.0E-5);
    ASSERT_LT(fabs(state.v.z-state2.v.z), 1.0E-5);
}

TEST_F(OrbitTest, OrbitElementsToStateVectorSpice1)
{
    double rpd = astro::RADPERDEG;

    // [1], example 4.7: 
    // Set the orbit elements directly    
    astro::OrbitElements oe;
    oe.h = 80000; // km²/s
    oe.e = 1.4; 
    oe.i = 30.0 * rpd;
    oe.omega = 40.0 * rpd;
    oe.w = 60.0 * rpd;
    oe.theta = 30 * rpd;
    oe.epoch = et;
    oe.mu = mu_earth;
    oe.rp = oe.h*oe.h / oe.mu * (1.0 / (1.0 + oe.e));
    oe.M0 = astro::OrbitElements::meanAnomalyFromTrueAnomaly(oe.theta, oe.e);
    print(oe);

    astro::State state0 = oe.toStateVectorOE();
	astro::State state1 = oe.toStateVectorSpice(et);
    
    print(state0);
    print(state1);

    astro::OrbitElements oe2 = astro::OrbitElements::fromStateVectorOE(state1, et, mu_earth);
    print(oe2);
    astro::OrbitElements oe3 = astro::OrbitElements::fromStateVectorSpice(state1, et, mu_earth);
    print(oe3);
    //ASSERT_EQ(1,0);

}
TEST_F(OrbitTest, OrbitElementsToStateVectorSpice2)
{
    // [1], Example 4.3:
    astro::State   state;
    state.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
    state.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]

    astro::OrbitElements oe = astro::OrbitElements::fromStateVectorSpice(state, et, mu_earth);
	print(oe);

	astro::State state1 = oe.toStateVectorSpice(et);
	std::cout << "Original" << std::endl;
	print(state);
	std::cout << "Recovered" << std::endl;
    print(state1);






}
