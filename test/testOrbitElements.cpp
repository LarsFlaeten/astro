// Reference from some of the tests:
// [1]  Orbital Mechanics for Engineering Students, 2nd Edition, Howard D. Curtis


#include "../astro/OrbitElements.cpp"
#include "../astro/State.h"
#include <gtest/gtest.h>

#include <cmath>

// Prettyprinters for vectors etc
#include "OrkExt.h"

using ork::vec3d;
using ork::mat3d;
using ork::mat4d;

using namespace astro;

class OrbitElementsTest : public ::testing::Test {

protected:
    OrbitElementsTest();

    virtual ~OrbitElementsTest();

    // Code here will be called immediately after the constructor (right
    // before each test).
    virtual void SetUp();

    // Code here will be called immediately after each test (right
    // before the destructor).
    virtual void TearDown();


    EphemerisTime et;
    double mu_earth;
};



OrbitElementsTest::OrbitElementsTest()
   : et(), mu_earth(398600.0)
{
    
}

OrbitElementsTest::~OrbitElementsTest()
{

}

void OrbitElementsTest::SetUp()
{
}

void OrbitElementsTest::TearDown()
{
}

TEST_F(OrbitElementsTest, OrbitalElementsFromState)
{

    // [1], Example 4.3:
    astro::PosState   state;
    state.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
    state.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]

    astro::OrbitElements oe = astro::OrbitElements::fromStateVectorOE(state, et, mu_earth);
    // Compare with [1] (Those are quite rounded numbers, so we cant be too accurate)
    ASSERT_LT(fabs(oe.h- 58310), 2.0);
    ASSERT_LT(fabs(oe.i * astro::DEGPERRAD - 153.2), 1.0E-1);
    ASSERT_LT(fabs(oe.omega * astro::DEGPERRAD - 255.3), 1.0E-1);
    ASSERT_LT(fabs(oe.e - 0.1712), 1.0E-4);
    ASSERT_LT(fabs(oe.w * astro::DEGPERRAD - 20.07), 1.0E-2);
    double theta = astro::OrbitElements::trueAnomalyFromMeanAnomaly(oe.M0, oe.e);

    ASSERT_LT(fabs(theta * astro::DEGPERRAD - 28.45), 1.0E-2);

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

TEST_F(OrbitElementsTest, OrbitalElementsFromState2)
{

    // [1], Example 4.3:
    astro::PosState   state;
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
    ASSERT_LT(fabs(oe.M0-oe2.M0), 0.0001);
    ASSERT_LT(fabs(oe.rp-oe2.rp), 0.0001);
    
    //print(oe);   
    //print(oe2);
}

TEST_F(OrbitElementsTest, OrbitalElementsFromStateCornerCases)
{
	double dpr = astro::DEGPERRAD;
    double mu = mu_earth;

    // [1], Example 4.3:
    astro::PosState   state;
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
    ASSERT_LT(fabs(oe.M0-oe2.M0), 0.0001);
    ASSERT_LT(fabs(oe.rp-oe2.rp), 0.0001);
    ASSERT_LT(fabs(oe.ap-oe2.ap), 0.0001);
    ASSERT_LT(fabs(oe.a-oe2.a), 0.0001);
    ASSERT_LT(fabs(oe.T-oe2.T), 0.0001);
    ASSERT_LT(fabs(oe.n-oe2.n), 0.0001);
      
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
    ASSERT_LT(fabs(oe.M0-oe2.M0), 0.0001);
    ASSERT_LT(fabs(oe.rp-oe2.rp), 0.0001);
    ASSERT_LT(fabs(oe.ap-oe2.ap), 0.0001);
    ASSERT_LT(fabs(oe.a-oe2.a), 0.0001);
    ASSERT_LT(fabs(oe.T-oe2.T), 0.0001);
    ASSERT_LT(fabs(oe.n-oe2.n), 0.0001);
 
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
    ASSERT_LT(fabs(oe.M0-oe2.M0), 0.0001);
    ASSERT_LT(fabs(oe.rp-oe2.rp), 0.0001);
    ASSERT_LT(fabs(oe.ap-oe2.ap), 0.0001);
    ASSERT_LT(fabs(oe.a-oe2.a), 0.0001);
    ASSERT_LT(fabs(oe.T-oe2.T), 0.0001);
    ASSERT_LT(fabs(oe.n-oe2.n), 0.0001);
 
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
    ASSERT_LT(fabs(oe.M0-oe2.M0), 0.0001);
    ASSERT_LT(fabs(oe.rp-oe2.rp), 0.0001);
    ASSERT_LT(fabs(oe.ap-oe2.ap), 0.0001);
    ASSERT_LT(fabs(oe.a-oe2.a), 0.0001);
    ASSERT_LT(fabs(oe.T-oe2.T), 0.0001);
    ASSERT_LT(fabs(oe.n-oe2.n), 0.0001);
 

    ASSERT_LT(fabs(oe.h-80000), 30);
    ASSERT_LT(fabs(oe.i*dpr-30), 0.01);
    ASSERT_LT(fabs(oe.omega*dpr-40), 0.022);
    ASSERT_LT(fabs(oe.e-1.4), 1.0E-2);
    ASSERT_LT(fabs(oe.w*dpr-60), 0.02);






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
    ASSERT_LT(fabs(oe.M0-oe2.M0), 0.0001);
    ASSERT_LT(fabs(oe.rp-oe2.rp), 0.0001);
    ASSERT_LT(fabs(oe.ap-oe2.ap), 0.0001);
    ASSERT_LT(fabs(oe.a-oe2.a), 0.0001);
    ASSERT_LT(fabs(oe.T-oe2.T), 0.0001);
    ASSERT_LT(fabs(oe.n-oe2.n), 0.0001);
 
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
    ASSERT_LT(fabs(oe.M0-oe2.M0), 0.0001);
    ASSERT_LT(fabs(oe.rp-oe2.rp), 0.0001);
	ASSERT_LT(fabs(oe.ap-oe2.ap), 0.0001);
    ASSERT_LT(fabs(oe.a-oe2.a), 0.0001);
    ASSERT_LT(fabs(oe.T-oe2.T), 0.0001);
    ASSERT_LT(fabs(oe.n-oe2.n), 0.0001);
 
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
    ASSERT_LT(fabs(oe.M0-oe2.M0), 0.0001);
    ASSERT_LT(fabs(oe.rp-oe2.rp), 0.0001);
    ASSERT_LT(fabs(oe.ap-oe2.ap), 0.0001);
    ASSERT_LT(fabs(oe.a-oe2.a), 0.0001);
    ASSERT_LT(fabs(oe.T-oe2.T), 0.0001);
    ASSERT_LT(fabs(oe.n-oe2.n), 0.0001);
 

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
    ASSERT_LT(fabs(oe.M0-oe2.M0), 0.0001);
    ASSERT_LT(fabs(oe.rp-oe2.rp), 0.0001);
    ASSERT_LT(fabs(oe.ap-oe2.ap), 0.0001);
    ASSERT_LT(fabs(oe.a-oe2.a), 0.0001);
    ASSERT_LT(fabs(oe.T-oe2.T), 0.0001);
    ASSERT_LT(fabs(oe.n-oe2.n), 0.0001);
 
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
    ASSERT_LT(fabs(oe.M0-oe2.M0), 0.0001);
    ASSERT_LT(fabs(oe.rp-oe2.rp), 0.0001);
    ASSERT_LT(fabs(oe.ap-oe2.ap), 0.0001);
    ASSERT_LT(fabs(oe.a-oe2.a), 0.0001);
    ASSERT_LT(fabs(oe.T-oe2.T), 0.0001);
    ASSERT_LT(fabs(oe.n-oe2.n), 0.0001);
 
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
    ASSERT_LT(fabs(oe.M0-oe2.M0), 0.0001);
    ASSERT_LT(fabs(oe.rp-oe2.rp), 0.0001);
    ASSERT_LT(fabs(oe.ap-oe2.ap), 0.0001);
    ASSERT_LT(fabs(oe.a-oe2.a), 0.0001);
    ASSERT_LT(fabs(oe.T-oe2.T), 0.0001);
    ASSERT_LT(fabs(oe.n-oe2.n), 0.0001);
 




    
}

#if 0
TEST_F(OrbitElementsTest, OrbitalElementsFromStateOpt)
{
    double mu_earth = 398600.0;

    // [1], Example 4.3:
    astro::PosState   state;
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
TEST_F(OrbitElementsTest, OrbitalElementsOEBenchmarkTime)
{
    // [1], Example 4.3:
    astro::PosState   state;
    state.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
    state.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]

    astro::OrbitElements oe;
    for(auto i = 0; i < 1000000; ++i)
    {
        oe = astro::OrbitElements::fromStateVectorOE(state, et, mu_earth);
    }

}

TEST_F(OrbitElementsTest, OrbitalElementsSpiceBenchmarkTime)
{

    // [1], Example 4.3:
    astro::PosState   state;
    state.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
    state.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]

    astro::OrbitElements oe;
    for(auto i = 0; i < 1000000; ++i)
    {
        oe = astro::OrbitElements::fromStateVectorSpice(state, et, mu_earth);
    }

}

#if 0
TEST_F(OrbitElementsTest, OrbitalElementsOptimizedTime)
{
    double mu_earth = 398600.0;

    // [1], Example 4.3:
    astro::PosState   state;
    state.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
    state.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]

    astro::OrbitElements oe;
    for(auto i = 0; i < 10000000; ++i)
    {
        oe = astro::OrbitElements::fromstateVectorOEOpt(state, mu_earth);
    }

}
#endif

TEST_F(OrbitElementsTest, OrkRotationTest)
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

TEST_F(OrbitElementsTest, OrbitElementsToStateVector)
{
    double rpd = astro::RADPERDEG;

    // [1], example 4.7:
    astro::OrbitElements oe;
    oe.h = 80000; // km²/s
    oe.e = 1.4;
    oe.i = 30.0 * rpd;
    oe.omega = 40.0 * rpd;
    oe.w = 60.0 * rpd;
    oe.M0 = astro::OrbitElements::meanAnomalyFromTrueAnomaly(30 * rpd, oe.e);
    oe.epoch = et;
    oe.mu = mu_earth;
    oe.computeDerivedQuantities();

    astro::PosState state = oe.toStateVectorOE(et);


    astro::OrbitElements oe2 = astro::OrbitElements::fromStateVectorOE(state, et, mu_earth);


    ASSERT_LT(fabs(oe2.h-oe.h), 1.0E-5);
    ASSERT_LT(fabs(oe2.i-oe.i), 1.0E-5);
    ASSERT_LT(fabs(oe2.omega-oe.omega), 1.0E-5);
    ASSERT_LT(fabs(oe2.w-oe.w), 1.0E-5);
    ASSERT_LT(fabs(oe2.e-oe.e), 1.0E-5);
    ASSERT_LT(fabs(oe.rp-oe2.rp), 0.0001);
    ASSERT_LT(fabs(oe.ap-oe2.ap), 0.0001);
    ASSERT_LT(fabs(oe.a-oe2.a), 0.0001);
    ASSERT_LT(fabs(oe.T-oe2.T), 0.0001);
    ASSERT_LT(fabs(oe.n-oe2.n), 0.0001);
 

}

TEST_F(OrbitElementsTest, OrbitElementsToStateVector2)
{
    // [1], Example 4.3:
    astro::PosState   state;
    state.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
    state.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]

    astro::OrbitElements oe = astro::OrbitElements::fromStateVectorOE(state, et, mu_earth);

    astro::PosState state2 = oe.toStateVectorOE(et);
    //print(state);
    //print(state2);
	ASSERT_LT(fabs(state.r.x-state2.r.x), 1.0E-5);
    ASSERT_LT(fabs(state.r.y-state2.r.y), 1.0E-5);
    ASSERT_LT(fabs(state.r.z-state2.r.z), 1.0E-5);
    ASSERT_LT(fabs(state.v.x-state2.v.x), 1.0E-5);
    ASSERT_LT(fabs(state.v.y-state2.v.y), 1.0E-5);
    ASSERT_LT(fabs(state.v.z-state2.v.z), 1.0E-5);
}

TEST_F(OrbitElementsTest, OrbitElementsToStateVectorSpice1)
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
    oe.epoch = et;
    oe.mu = mu_earth;
    oe.M0 = astro::OrbitElements::meanAnomalyFromTrueAnomaly(30.0*rpd, oe.e);
    oe.computeDerivedQuantities();

    //print(oe);

    astro::PosState state0 = oe.toStateVectorOE(et);
	astro::PosState state1 = oe.toStateVectorSpice(et);
    
    //print(state0);
    //print(state1);

    astro::OrbitElements oe2 = astro::OrbitElements::fromStateVectorOE(state0, et, mu_earth);
    //print(oe2);
    astro::OrbitElements oe3 = astro::OrbitElements::fromStateVectorSpice(state1, et, mu_earth);
    //print(oe3);
    //ASSERT_EQ(1,0);

    ASSERT_LT(fabs(oe2.h-oe.h), 1.0E-5);
    ASSERT_LT(fabs(oe2.i-oe.i), 1.0E-5);
    ASSERT_LT(fabs(oe2.omega-oe.omega), 1.0E-5);
    ASSERT_LT(fabs(oe2.w-oe.w), 1.0E-5);
    ASSERT_LT(fabs(oe2.e-oe.e), 1.0E-5);
	ASSERT_LT(fabs(oe2.epoch.getETValue()-oe.epoch.getETValue()), 1.0E-5);
	ASSERT_LT(fabs(oe2.mu-oe.mu), 1.0E-5);
	ASSERT_LT(fabs(oe2.rp-oe.rp), 1.0E-5);
	ASSERT_LT(fabs(oe2.M0-oe.M0), 1.0E-5);
    ASSERT_LT(fabs(oe.ap-oe2.ap), 0.0001);
    ASSERT_LT(fabs(oe.a-oe2.a), 0.0001);
    ASSERT_LT(fabs(oe.T-oe2.T), 0.0001);
    ASSERT_LT(fabs(oe.n-oe2.n), 0.0001);
 
    ASSERT_LT(fabs(oe3.h-oe.h), 1.0E-5);
    ASSERT_LT(fabs(oe3.i-oe.i), 1.0E-5);
    ASSERT_LT(fabs(oe3.omega-oe.omega), 1.0E-5);
    ASSERT_LT(fabs(oe3.w-oe.w), 1.0E-5);
    ASSERT_LT(fabs(oe3.e-oe.e), 1.0E-5);
    ASSERT_LT(fabs(oe3.epoch.getETValue()-oe.epoch.getETValue()), 1.0E-5);
    ASSERT_LT(fabs(oe3.mu-oe.mu), 1.0E-5);
    ASSERT_LT(fabs(oe3.rp-oe.rp), 1.0E-5);
    ASSERT_LT(fabs(oe3.M0-oe.M0), 1.0E-5);
    ASSERT_LT(fabs(oe.ap-oe2.ap), 0.0001);
    ASSERT_LT(fabs(oe.a-oe2.a), 0.0001);
    ASSERT_LT(fabs(oe.T-oe2.T), 0.0001);
    ASSERT_LT(fabs(oe.n-oe2.n), 0.0001);
 



}
TEST_F(OrbitElementsTest, OrbitElementsToStateVectorSpice2)
{
    // [1], Example 4.3:
    astro::PosState   state;
    state.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
    state.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]

    astro::OrbitElements oe = astro::OrbitElements::fromStateVectorSpice(state, et, mu_earth);
	//print(oe);

	astro::PosState state1 = oe.toStateVectorSpice(et);
	//std::cout << "Original" << std::endl;
	//print(state);
	//std::cout << "Recovered" << std::endl;
    //print(state1);

    ASSERT_LT(fabs(state.r.x-state1.r.x), 1.0E-5);
    ASSERT_LT(fabs(state.r.y-state1.r.y), 1.0E-5);
    ASSERT_LT(fabs(state.r.z-state1.r.z), 1.0E-5);
    ASSERT_LT(fabs(state.v.x-state1.v.x), 1.0E-5);
    ASSERT_LT(fabs(state.v.y-state1.v.y), 1.0E-5);
    ASSERT_LT(fabs(state.v.z-state1.v.z), 1.0E-5);

}




TEST_F(OrbitElementsTest, AnomalyConversions)
{

    

	    
    for(double v = 0.0; v < astro::TWOPI; v+= astro::PIHALF/10.0)
    {
        for(double e = 0.0; e < 0.999; e+= 0.1)
		{
            double M = astro::OrbitElements::meanAnomalyFromTrueAnomaly(v, e);
            double v2 = astro::OrbitElements::trueAnomalyFromMeanAnomaly(M, e);
            ASSERT_LT(fabs(v-v2), astro::KEPLER_TOLERANCE);
        }

        // Upper level of eccentricity for Kepler2 to be accurate
        // is around 0.9995 for all theta = [0..2pi]
        double e = 0.99979;
        double M = astro::OrbitElements::meanAnomalyFromTrueAnomaly(v, e);
        double v2 = astro::OrbitElements::trueAnomalyFromMeanAnomaly(M, e);
        //std::cout << v << ", " << e << std::endl;
        //std::cout << fabs(v-v2) << ", " << v << std::endl;

        // We allow a little bit higher tolerance on true anomaly
        // compared to eccentric anomaly, since it is a derived value
        ASSERT_LT(fabs(v-v2), astro::KEPLER_TOLERANCE*5.0);
    }


}

TEST_F(OrbitElementsTest, KeplerCompTest1)
{
    double rpd = astro::RADPERDEG;
    double M = 30.0*rpd;
    double e = 0.6;
    std::pair<double,int> res1 = astro::OrbitElements::Kepler1(M, e);
    //std::cout << "Kepler1" << std::endl;
    //std::cout << "  Mean anomaly: " << M << std::endl;
    //std::cout << "  Eccentricity: " << e << std::endl;
    //std::cout << "  Ecc. Anomaly: " << res1.first << std::endl;
    //std::cout << "  Iterations:   " << res1.second << std::endl;
	// Check on anomaly:
    double E1 = res1.first;
    double M12 = E1 - e*sin(E1);
    ASSERT_LT(fabs(M - M12), astro::KEPLER_TOLERANCE); 

    std::pair<double,int> res2 = astro::OrbitElements::Kepler2(M, e);
    //std::cout << "Kepler2" << std::endl;
    //std::cout << "  Mean anomaly: " << M << std::endl;
    //std::cout << "  Eccentricity: " << e << std::endl;
    //std::cout << "  Ecc. Anomaly: " << res2.first << std::endl;
    //std::cout << "  Iterations:   " << res2.second << std::endl;
	// Check on anomaly:
    double E2 = res1.first;
    double M22 = E2 - e*sin(E2);
    ASSERT_LT(fabs(M - M22), astro::KEPLER_TOLERANCE); 

    // Check the two independent methods
    ASSERT_LT(fabs(res1.first - res2.first), astro::KEPLER_TOLERANCE); 
}

// Results from running below with release build (-O2)
//[ RUN      ] OrbitElementsTest.Kepler1BenchMark
//[       OK ] OrbitElementsTest.Kepler1BenchMark (1660 ms)
//[ RUN      ] OrbitElementsTest.Kepler2BenchMark
//[       OK ] OrbitElementsTest.Kepler2BenchMark (226 ms)

TEST_F(OrbitElementsTest, Kepler1BenchMark)
{
    double rpd = astro::RADPERDEG;
    double M = 30.0*rpd;
    double e = 0.6;

    for(int i = 0; i < 1000000; i++)
        std::pair<double,int> res = astro::OrbitElements::Kepler1(M,e);

}

TEST_F(OrbitElementsTest, Kepler2BenchMark)
{
    double rpd = astro::RADPERDEG;
    double M = 30.0*rpd;
    double e = 0.6;

    for(int i = 0; i < 1000000; i++)
        std::pair<double,int> res = astro::OrbitElements::Kepler2(M,e);

}

TEST_F(OrbitElementsTest, Kepler2HypBenchMark)
{
    double rpd = astro::RADPERDEG;
    double M = 30.0*rpd;
    double e = 1.4;

    for(int i = 0; i < 1000000; i++)
        std::pair<double,int> res = astro::OrbitElements::Kepler2(M,e);

}



TEST_F(OrbitElementsTest, Kepler2ConvergenceTest)
{
    double rpd = astro::RADPERDEG;
	double M, M2;
    for(double m = 0.0; m < 360.0; m = m + 10)
		for(double e = 0.0; e < 1.0; e += 0.1)
		{	

            if ( e>= 0.9998)
                e = 0.9997;
			M = m*rpd;
        	std::pair<double,int> res = astro::OrbitElements::Kepler2(M,e);
			
            // Check on anomaly:
            double E = res.first;
            M2 = E - e*sin(E);
            ASSERT_LT(fabs(M - M2), astro::KEPLER_TOLERANCE); 
            // Uncomment this line for report on convergence:
			// std::cout << "it: " << res.second << ", M: " << M << ", e: " << e << ", E: " << res.first << std::endl;
		}
}

TEST_F(OrbitElementsTest, Kepler2HypTest1)
{
   	double rpd = astro::RADPERDEG;
    double M = -30.0*rpd;
    double e = 1.4;

    std::pair<double,int> res;
	
    ASSERT_THROW(res = astro::OrbitElements::Kepler1(M, e), astro::AstroException);

    res = astro::OrbitElements::Kepler2(M, e);
    double H = res.first;
    double M2 = e*sinh(H)-H;
	//std::cout << "M: " << M << ", M2: " << M2 << ", e: " << e << ", H: " << res.first << ", it: " << res.second << std::endl;
    ASSERT_LT(fabs(M-M2), astro::KEPLER_TOLERANCE);

	e = 1.1;
	res = astro::OrbitElements::Kepler2(M, e);
    H = res.first;
    M2 = e*sinh(H)-H;
	//std::cout << "M: " << M << ", M2: " << M2 << ", e: " << e << ", H: " << res.first << ", it: " << res.second << std::endl;
    ASSERT_LT(fabs(M-M2), astro::KEPLER_TOLERANCE);


    e = 1.05;
    res = astro::OrbitElements::Kepler2(M, e);
    H = res.first;
    M2 = e*sinh(H)-H;
	//std::cout << "M: " << M << ", M2: " << M2 << ", e: " << e << ", H: " << res.first << ", it: " << res.second << std::endl;
    ASSERT_LT(fabs(M-M2), astro::KEPLER_TOLERANCE);


    e = 2.1;
    res = astro::OrbitElements::Kepler2(M, e);
    H = res.first;
    M2 = e*sinh(H)-H;
	//std::cout << "M: " << M << ", M2: " << M2 << ", e: " << e << ", H: " << res.first << ", it: " << res.second << std::endl;
    ASSERT_LT(fabs(M-M2), astro::KEPLER_TOLERANCE);


	e = 3.1;
    res = astro::OrbitElements::Kepler2(M, e);
    H = res.first;
    M2 = e*sinh(H)-H;
	//std::cout << "M: " << M << ", M2: " << M2 << ", e: " << e << ", H: " << res.first << ", it: " << res.second << std::endl;
    ASSERT_LT(fabs(M-M2), astro::KEPLER_TOLERANCE);


    e = 1.4;
    M = 0.0*rpd;
    res = astro::OrbitElements::Kepler2(M, e);
    H = res.first;
    M2 = e*sinh(H)-H;
	//std::cout << "M: " << M << ", M2: " << M2 << ", e: " << e << ", H: " << res.first << ", it: " << res.second << std::endl;
    ASSERT_LT(fabs(M-M2), astro::KEPLER_TOLERANCE);

    M = -10.0*rpd;
    res = astro::OrbitElements::Kepler2(M, e);
    H = res.first;
    M2 = e*sinh(H)-H;
	//std::cout << "M: " << M << ", M2: " << M2 << ", e: " << e << ", H: " << res.first << ", it: " << res.second << std::endl;
    ASSERT_LT(fabs(M-M2), astro::KEPLER_TOLERANCE);


    M = -20.0*rpd;
    res = astro::OrbitElements::Kepler2(M, e);
    H = res.first;
    M2 = e*sinh(H)-H;
	//std::cout << "M: " << M << ", M2: " << M2 << ", e: " << e << ", H: " << res.first << ", it: " << res.second << std::endl;
    ASSERT_LT(fabs(M-M2), astro::KEPLER_TOLERANCE);


    M = -30.0*rpd;
    res = astro::OrbitElements::Kepler2(M, e);
    H = res.first;
    M2 = e*sinh(H)-H;
	//std::cout << "M: " << M << ", M2: " << M2 << ", e: " << e << ", H: " << res.first << ", it: " << res.second << std::endl;
    ASSERT_LT(fabs(M-M2), astro::KEPLER_TOLERANCE);


    M = -40.0*rpd;
    res = astro::OrbitElements::Kepler2(M, e);
    H = res.first;
    M2 = e*sinh(H)-H;
	//std::cout << "M: " << M << ", M2: " << M2 << ", e: " << e << ", H: " << res.first << ", it: " << res.second << std::endl;
    ASSERT_LT(fabs(M-M2), astro::KEPLER_TOLERANCE);


    M = -50.0*rpd;
    res = astro::OrbitElements::Kepler2(M, e);
    H = res.first;
    M2 = e*sinh(H)-H;
	//std::cout << "M: " << M << ", M2: " << M2 << ", e: " << e << ", H: " << res.first << ", it: " << res.second << std::endl;
    ASSERT_LT(fabs(M-M2), astro::KEPLER_TOLERANCE);


    M = -60.0*rpd;
    res = astro::OrbitElements::Kepler2(M, e);
    H = res.first;
    M2 = e*sinh(H)-H;
	//std::cout << "M: " << M << ", M2: " << M2 << ", e: " << e << ", H: " << res.first << ", it: " << res.second << std::endl;
    ASSERT_LT(fabs(M-M2), astro::KEPLER_TOLERANCE);


    M = -70.0*rpd;
    res = astro::OrbitElements::Kepler2(M, e);
    H = res.first;
    M2 = e*sinh(H)-H;
	//std::cout << "M: " << M << ", M2: " << M2 << ", e: " << e << ", H: " << res.first << ", it: " << res.second << std::endl;
    ASSERT_LT(fabs(M-M2), astro::KEPLER_TOLERANCE);

    M = -80.0*rpd;
    res = astro::OrbitElements::Kepler2(M, e);
    H = res.first;
    M2 = e*sinh(H)-H;
	//std::cout << "M: " << M << ", M2: " << M2 << ", e: " << e << ", H: " << res.first << ", it: " << res.second << std::endl;
    ASSERT_LT(fabs(M-M2), astro::KEPLER_TOLERANCE);

    M = -90.0*rpd;
    res = astro::OrbitElements::Kepler2(M, e);
    H = res.first;
    M2 = e*sinh(H)-H;
	//std::cout << "M: " << M << ", M2: " << M2 << ", e: " << e << ", H: " << res.first << ", it: " << res.second << std::endl;
    ASSERT_LT(fabs(M-M2), astro::KEPLER_TOLERANCE);

    e = 1.05;
    M = -90.0*rpd;
    res = astro::OrbitElements::Kepler2(M, e);
    H = res.first;
    M2 = e*sinh(H)-H;
	//std::cout << "M: " << M << ", M2: " << M2 << ", e: " << e << ", H: " << res.first << ", it: " << res.second << std::endl;
    ASSERT_LT(fabs(M-M2), astro::KEPLER_TOLERANCE);




 
}

TEST_F(OrbitElementsTest, Kepler2HypConvergenceTest)
{
    double rpd = astro::RADPERDEG;
	double M, M2;
    for(double m = -90; m < 90.0; m = m + 10)
		for(double e = 1.05; e < 3.0; e += 0.1)
		{	

			M = m*rpd;
        	std::pair<double,int> res = astro::OrbitElements::Kepler2(M,e);
			
            // Check on anomaly:
            double H = res.first;
            M2 = e*sinh(H)-H;
            // Uncomment this line for report on convergence:
			// std::cout << "it: " << res.second << ", M: " << M << ", Mcheck: " << M2 << ", e: " << e << ", H: " << res.first << std::endl;
	        ASSERT_LT(fabs(M - M2), astro::KEPLER_TOLERANCE); 
 	}
}

// The formal test of OE-methods against SPICE before we do benchmarking
TEST_F(OrbitElementsTest, OrbitElementsToStateVectorComp)
{
    // [1], Example 4.3:
    astro::PosState   state;
    state.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
    state.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]

    astro::OrbitElements oe1 = astro::OrbitElements::fromStateVectorOE(state, et, mu_earth);

    // [1], example 4.7: 
    // Set the orbit elements directly    
    double rpd = astro::RADPERDEG;

    astro::OrbitElements oe2;
    oe2.h = 80000; // km²/s
    oe2.e = 1.4; 
    oe2.i = 30.0 * rpd;
    oe2.omega = 40.0 * rpd;
    oe2.w = 60.0 * rpd;
    oe2.epoch = et;
    oe2.mu = mu_earth;
    oe2.M0 = astro::OrbitElements::meanAnomalyFromTrueAnomaly(30.0*rpd, oe2.e);
    oe2.computeDerivedQuantities();



    // Get state 1000s later:
    EphemerisTime et2(1000.0);
    astro::PosState state1 = oe1.toStateVectorOE(et2);
    astro::PosState state2 = oe1.toStateVectorSpice(et2);
    //print(state1);
    //print(state2);
    ASSERT_LT(fabs(state2.r.x-state1.r.x), 1.0E-5);
    ASSERT_LT(fabs(state2.r.y-state1.r.y), 1.0E-5);
    ASSERT_LT(fabs(state2.r.z-state1.r.z), 1.0E-5);
    ASSERT_LT(fabs(state2.v.x-state1.v.x), 1.0E-5);
    ASSERT_LT(fabs(state2.v.y-state1.v.y), 1.0E-5);
    ASSERT_LT(fabs(state2.v.z-state1.v.z), 1.0E-5);

    astro::PosState state3 = oe2.toStateVectorOE(et2);
    astro::PosState state4 = oe2.toStateVectorSpice(et2);
    ASSERT_LT(fabs(state3.r.x-state4.r.x), 1.0E-5);
    ASSERT_LT(fabs(state3.r.y-state4.r.y), 1.0E-5);
    ASSERT_LT(fabs(state3.r.z-state4.r.z), 1.0E-5);
    ASSERT_LT(fabs(state3.v.x-state4.v.x), 1.0E-5);
    ASSERT_LT(fabs(state3.v.y-state4.v.y), 1.0E-5);
    ASSERT_LT(fabs(state3.v.z-state4.v.z), 1.0E-5);
}



// Benchmark testing run on release (-O2), shows that OE is faster
// [==========] Running 4 tests from 1 test case.
// [----------] Global test environment set-up.
// [----------] 4 tests from OrbitElementsTest
// [ RUN      ] OrbitElementsTest.OrbitElementsToStateVectorBenchmarkEllipticOE
// [       OK ] OrbitElementsTest.OrbitElementsToStateVectorBenchmarkEllipticOE (546 ms)
// [ RUN      ] OrbitElementsTest.OrbitElementsToStateVectorBenchmarkEllipticSpice
// [       OK ] OrbitElementsTest.OrbitElementsToStateVectorBenchmarkEllipticSpice (3751 ms)
// [ RUN      ] OrbitElementsTest.OrbitElementsToStateVectorBenchmarkHypOE
// [       OK ] OrbitElementsTest.OrbitElementsToStateVectorBenchmarkHypOE (1195 ms)
// [ RUN      ] OrbitElementsTest.OrbitElementsToStateVectorBenchmarkHypSpice
// [       OK ] OrbitElementsTest.OrbitElementsToStateVectorBenchmarkHypSpice (3076 ms)
// [----------] 4 tests from OrbitElementsTest (8568 ms total)
//
// [----------] Global test environment tear-down
// [==========] 4 tests from 1 test case ran. (8568 ms total)
// [  PASSED  ] 4 tests.

TEST_F(OrbitElementsTest, OrbitElementsToStateVectorBenchmarkEllipticOE)
{
    // [1], Example 4.3:
    astro::PosState   state;
    state.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
    state.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]

    astro::OrbitElements oe1 = astro::OrbitElements::fromStateVectorOE(state, et, mu_earth);


    // Get state 1000s later:
    EphemerisTime et2(1000.0);

    for(int i = 0; i < 1000000; ++i)
        astro::PosState state1 = oe1.toStateVectorOE(et2);
}

TEST_F(OrbitElementsTest, OrbitElementsToStateVectorBenchmarkEllipticSpice)
{
    // [1], Example 4.3:
    astro::PosState   state;
    state.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
    state.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]

    astro::OrbitElements oe1 = astro::OrbitElements::fromStateVectorOE(state, et, mu_earth);



    // Get state 1000s later:
    EphemerisTime et2(1000.0);

    for(int i = 0; i < 1000000; ++i)
        astro::PosState state1 = oe1.toStateVectorSpice(et2);
}

TEST_F(OrbitElementsTest, OrbitElementsToStateVectorBenchmarkHypOE)
{


    // [1], example 4.7: 
    // Set the orbit elements directly    
    double rpd = astro::RADPERDEG;

    astro::OrbitElements oe2;
    oe2.h = 80000; // km²/s
    oe2.e = 1.4; 
    oe2.i = 30.0 * rpd;
    oe2.omega = 40.0 * rpd;
    oe2.w = 60.0 * rpd;
    oe2.epoch = et;
    oe2.mu = mu_earth;
    oe2.M0 = astro::OrbitElements::meanAnomalyFromTrueAnomaly(30.0*rpd, oe2.e);
    oe2.computeDerivedQuantities();



    // Get state 1000s later:
    EphemerisTime et2(1000.0);

    for(int i = 0; i < 1000000; ++i)
        astro::PosState state1 = oe2.toStateVectorOE(et2);
}

TEST_F(OrbitElementsTest, OrbitElementsToStateVectorBenchmarkHypSpice)
{


    // [1], example 4.7: 
    // Set the orbit elements directly    
    double rpd = astro::RADPERDEG;

    astro::OrbitElements oe2;
    oe2.h = 80000; // km²/s
    oe2.e = 1.4; 
    oe2.i = 30.0 * rpd;
    oe2.omega = 40.0 * rpd;
    oe2.w = 60.0 * rpd;
    oe2.epoch = et;
    oe2.mu = mu_earth;
    oe2.M0 = astro::OrbitElements::meanAnomalyFromTrueAnomaly(30.0*rpd, oe2.e);
    oe2.computeDerivedQuantities();



    // Get state 1000s later:
    EphemerisTime et2(1000.0);

    for(int i = 0; i < 1000000; ++i)
        astro::PosState state1 = oe2.toStateVectorSpice(et2);
}


TEST_F(OrbitElementsTest, OrbitPropagationTest1)
{
    // [1], Example 4.3:
    astro::PosState   state;
    state.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
    state.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]

    astro::OrbitElements oe1 = astro::OrbitElements::fromStateVector(state, et, mu_earth);


    // Get state one period later, should give same state
    double T = oe1.T;

    EphemerisTime et2(T);

    astro::PosState state1 = oe1.toStateVector(et2);
    ASSERT_LT(fabs(state.r.x-state1.r.x), 1.0E-5);
    ASSERT_LT(fabs(state.r.y-state1.r.y), 1.0E-5);
    ASSERT_LT(fabs(state.r.z-state1.r.z), 1.0E-5);
    ASSERT_LT(fabs(state.v.x-state1.v.x), 1.0E-5);
    ASSERT_LT(fabs(state.v.y-state1.v.y), 1.0E-5);
    ASSERT_LT(fabs(state.v.z-state1.v.z), 1.0E-5);

    // Get state 150 periods later, should give same state
    EphemerisTime et4(150*T);
    state1 = oe1.toStateVector(et4);
    ASSERT_LT(fabs(state.r.x-state1.r.x), 1.0E-5);
    ASSERT_LT(fabs(state.r.y-state1.r.y), 1.0E-5);
    ASSERT_LT(fabs(state.r.z-state1.r.z), 1.0E-5);
    ASSERT_LT(fabs(state.v.x-state1.v.x), 1.0E-5);
    ASSERT_LT(fabs(state.v.y-state1.v.y), 1.0E-5);
    ASSERT_LT(fabs(state.v.z-state1.v.z), 1.0E-5);

    // Get state 40 periods before, should give same state
    EphemerisTime et5(-40*T);
    state1 = oe1.toStateVector(et5);
    ASSERT_LT(fabs(state.r.x-state1.r.x), 1.0E-5);
    ASSERT_LT(fabs(state.r.y-state1.r.y), 1.0E-5);
    ASSERT_LT(fabs(state.r.z-state1.r.z), 1.0E-5);
    ASSERT_LT(fabs(state.v.x-state1.v.x), 1.0E-5);
    ASSERT_LT(fabs(state.v.y-state1.v.y), 1.0E-5);
    ASSERT_LT(fabs(state.v.z-state1.v.z), 1.0E-5);

    //print(oe1);

    // Get the time to periapsis and apoapsis:
    double k = oe1.T/astro::TWOPI;
    double e = oe1.e;

    double E_now = astro::OrbitElements::eccentricAnomalyFromMeanAnomaly(oe1.M0, oe1.e);
    double E_pe = 0.0;
    double E_ap = astro::PI;
    
    double t_now = k*(E_now - e*sin(E_now));
    double t_pe = k*(E_pe - e*sin(E_pe))-t_now;
    double t_ap = k*(E_ap - e*sin(E_ap))-t_now;
    
    t_pe = t_pe < t_now ? t_pe+T : t_pe;
    t_ap = t_ap < t_now ? t_ap+T : t_ap;



    //std::cout << "Period:      " << oe1.T << std::endl;
    //std::cout << "Time now:    " << et.getETValue() << std::endl;

    //std::cout << "Time to PE:  " << t_pe << std::endl;
    //std::cout << "Time to AP:  " << t_ap << std::endl;

    // Propagate the orbit to PE:
    EphemerisTime et_pe(t_pe);
    state1 = oe1.toStateVector(et_pe);
    double r_pe = state1.r.length();
    //std::cout << "Radius @ PE: " << r_pe << std::endl;
    ASSERT_LT(fabs(r_pe - oe1.rp), 1.0E-5);

    // Propagate the orbit to AP:
    EphemerisTime et_ap(t_ap);
    state1 = oe1.toStateVector(et_ap);
    double r_ap = state1.r.length();
    //std::cout << "Radius @ AP: " << r_ap << std::endl;
    ASSERT_LT(fabs(r_ap - oe1.ap), 1.0E-5);

    // Further tests om some propagations:

    // [1], example 3.5
    double rpd = astro::RADPERDEG;

   
    double v_pe = 15.0; // [km/s]
    r_pe = 300.0 + 6378.0;
    EphemerisTime et35(0.0);
    OrbitElements oe35;
    oe35.h = v_pe*r_pe;
    oe35.e = oe35.h*oe35.h/(mu_earth*r_pe)-1.0; // e = 2.7696
    oe35.i = 0.0;
    oe35.omega = 0.0;
    oe35.w = 0.0;
    oe35.M0 = 0.0;
    oe35.epoch = et35;
    oe35.mu = mu_earth;
    oe35.computeDerivedQuantities();
    //print(oe35);

    // Find radius and time when true anomaly is 100 degrees
    double v = 100*rpd;
    double Mh = astro::OrbitElements::meanAnomalyFromTrueAnomaly(v, oe35.e);
    
    double t100 = Mh/oe35.n;   
    ASSERT_LT(fabs(t100-4141), 1.0); // Time at v = 100 deg = 4141s

    PosState state35 = oe35.toStateVector(EphemerisTime(4141));
    double r100 = state35.r.length();
    //std::cout << r100 << std::endl;
    ASSERT_LT(fabs(r100-48497), 6); // Radius at 4141s = 48497 km
    // Note by going throug mean anomaly here, we get some rounding error
    // true anomaly (99.9989 vs 100). Together with rounding in ex3.5 on e
    // his will give a few kms off...

    // Find posiiton and speed 3 hours later
    double t3h = 4141 + 3*3600;
    state35 = oe35.toStateVector(EphemerisTime(t3h));
    ASSERT_LT(fabs(state35.r.length()-163180), 6); // Radius 3h later is 163180km
    ASSERT_LT(fabs(state35.v.length()-10.51), 1.0E-2); // velocity 3h later is 10.51km/s

    // [1], example 3.1
    EphemerisTime et31(0.0);
    OrbitElements oe31;
    oe31.h = 72472;
    oe31.e = 0.37255;
    oe31.i = 0.0;
    oe31.omega = 0.0;
    oe31.w = 0.0;
    oe31.M0 = 0.0;
    oe31.epoch = et35;
    oe31.mu = mu_earth;
    oe31.computeDerivedQuantities();
	print(oe31);

	// Calculate time to fly from PE to a true anomaly of 120 degrees:
	double M = astro::OrbitElements::meanAnomalyFromTrueAnomaly(120*rpd, oe31.e);
	double t31 = oe31.T*M/astro::TWOPI;
	 
	ASSERT_LT(fabs(t31-4077), 1.0E02);
	
	
}

TEST_F(OrbitElementsTest, HyperbolicAsymptote)
{
    astro::PosState s;
    s.r = vec3d(6378 + 400, 0.0, 0.0);
    //double v = 7.66895; // Circular
    double v = 10.84509; // roughly eccentric (to six decimals on e)
    s.v = vec3d(0, 0, v*1.1);
    

    astro::OrbitElements oe = astro::OrbitElements::fromStateVectorOE(s, et, mu_earth);
    //std::cout << oe << std::endl; 
    // Propagate the orbit for a long time:
    PosState s_inf = oe.toStateVector(et + 3.5*60*60*24);
    //std::cout << s << std::endl;
    //std::cout << s_inf << std::endl;
    
    // the state vector should now roughly be asymptotic.
    // Compare to hyperbolic excess velocity:
    // This one is a bit tricky, since velocity is slowly converging
    // and we will start to encounter problems in the Kepler solver
    // 3.5 days above seems to be eround the limit of cenvergence for
    // the above orbit..
    // TODO: We need to map this behaviour. 
    // Non-convergendce happens outside the SOI for this orbit (e = 1.42),
    // but this mght not be the case for other eccentricities.
    // Generally it is therfore unpredictable to base the simulator on
    // keplerian hyperbolic orbits for long time spans.. We might also
    // get into trouble with pacthed conics method. Will just have to wait and see.
    // Current maxiter is 100, we might have to relax on this..
    double v_inf = s_inf.v.length();
    ASSERT_LT(fabs(v_inf- astro::hyperbolicExcessVelocity(mu_earth, oe.a)), 0.06);
    
    // asymptotic angles are more easy
    double hyp_asym = astro::hyperbolicAsymptote(oe.e);
    vec3d v0 = -s.r.normalize();
    vec3d v1 = s_inf.v.normalize();
    vec3d v2 = s_inf.r.normalize();
    double hyp_asym1 = acos(v0.dotproduct(v1)); // Asymptote from velocity
    double hyp_asym2 = acos(v0.dotproduct(v2)); // Asymptote from position
  
    ASSERT_LT(fabs(hyp_asym - hyp_asym1), 0.0005); // 0.03 degrees..
    ASSERT_LT(fabs(hyp_asym - hyp_asym2), 0.02); // 0.03 degrees..


}
