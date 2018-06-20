#include "../astro/Util.h"
#include "../astro/Interpolate.cpp"
#include "../astro/ODE.h"
#include "../astro/OrbitElements.h"
#include "../astro/Propagator.h"
#include <gtest/gtest.h>

#include <cmath>


using namespace astro;

class InterpolateTest : public ::testing::Test {

protected:
    InterpolateTest();

    virtual ~InterpolateTest();

    // Code here will be called immediately after the constructor (right
    // before each test).
    virtual void SetUp();

    // Code here will be called immediately after each test (right
    // before the destructor).
    virtual void TearDown();
};



InterpolateTest::InterpolateTest()
{

}

InterpolateTest::~InterpolateTest()
{

}

void InterpolateTest::SetUp()
{
}

void InterpolateTest::TearDown()
{
}

// Test corner cases for different time spans:
TEST_F(InterpolateTest, CornerCases)
{
    PosState s1;
    s1.r = vec3d(0.0, 0.0, 0.0);
    s1.v = vec3d(1.0, 2.0, 3.0);
    EphemerisTime et1(1000);

    PosState s2;
    s2.r = vec3d(1.0, 2.0, 3.0);
    s2.v = vec3d(0.0, 0.0, 0.0);
    EphemerisTime et2(1001);

    PosState s3;
    astro::hermite(s1, et1, s2, et2, et1, s3);
    ASSERT_EQ(s3.r, vec3d(0.0, 0.0, 0.0));
    ASSERT_EQ(s3.v, vec3d(1.0, 2.0, 3.0));

    astro::hermite(s1, et1, s2, et2, et2, s3);
    ASSERT_EQ(s3.r, vec3d(1.0, 2.0, 3.0));
    ASSERT_EQ(s3.v, vec3d(0.0, 0.0, 0.0));
  

    et1 = EphemerisTime(0);
    et2 = EphemerisTime(1000);
    astro::hermite(s1, et1, s2, et2, et1, s3);
    ASSERT_EQ(s3.r, vec3d(0.0, 0.0, 0.0));
    ASSERT_EQ(s3.v, vec3d(1.0, 2.0, 3.0));

    astro::hermite(s1, et1, s2, et2, et2, s3);
    ASSERT_EQ(s3.r, vec3d(1.0, 2.0, 3.0));
    ASSERT_EQ(s3.v, vec3d(0.0, 0.0, 0.0));
    
    et1 = EphemerisTime(0);
    et2 = EphemerisTime(0.01);
    astro::hermite(s1, et1, s2, et2, et1, s3);
    ASSERT_EQ(s3.r, vec3d(0.0, 0.0, 0.0));
    ASSERT_EQ(s3.v, vec3d(1.0, 2.0, 3.0));

    astro::hermite(s1, et1, s2, et2, et2, s3);
    ASSERT_EQ(s3.r, vec3d(1.0, 2.0, 3.0));
    ASSERT_EQ(s3.v, vec3d(0.0, 0.0, 0.0));
 
}

TEST_F(InterpolateTest, Evolution)
{
    PosState s1;
    s1.r = vec3d(0.0, 0.0, 0.0);
    s1.v = vec3d(2.0, 4.0, 6.0);
    PosState s2;
    s2.r = vec3d(1000.0, 2000.0, 3000.0);
    s2.v = vec3d(0.0, 0.0, 0.0);
     
    PosState s3;
    astro::hermite(s1, EphemerisTime(0), s2, EphemerisTime(1000), EphemerisTime(500), s3);

    ASSERT_EQ(s3.v, vec3d(1, 2, 3));
    ASSERT_EQ(s3.r, vec3d(1000*3/4, 2000*3/4, 3000*3/4));
 
  
}

TEST_F(InterpolateTest, OrbitState1)
{
    astro::EphemerisTime et01(0);
    astro::EphemerisTime et02(2347589.0);
    double mu_earth(398600.0);
    astro::Attractor a = {vec3d::ZERO, mu_earth};
    astro::ODE  ode;
    ode.addAttractor(a);

    PosState state0;
    double v =  58311.7/7283.46;
    double r = 7283.46;
    state0.r = vec3d(r, 0.0 , 0.0);  //[km]
    state0.v = vec3d(0.0,1.2*v*sqrt(0.5), 0.5*v*sqrt(0.5));      //[km/s]

    OrbitElements oe01 = astro::OrbitElements::fromStateVector(state0, et01, mu_earth);

    OrbitElements oe02 = astro::OrbitElements::fromStateVector(state0, et02, mu_earth);


    astro::Propagator<astro::ODE, astro::RKF78> pr(ode);
    
    TimeDelta dt(1.0);
    EphemerisTime et = et01;    
    PosState state = state0;
    for(int i = 0; i < 6; i++)
    {
        auto res = pr.doStep(state, et, dt);
        dt = res.dt_next;
        et = res.et;
        state = res.s;
        //std::cout << dt.value << " - " << res.s.r << " - " << res.s.v << std::endl;
    }

    // Our interpolation interval from:
    auto res1 = pr.doStep(state, et, dt);

    // .. to:
    // Ok, this is important. By following RKF56/78s suggested dt_next
    // we get great intervals like 90s for RKF45 and 180s for RKF78
    // bu the interpolation is not that correct in the midpoint of these
    // RKF78s dt_next gave a few meters off in the mid point.
    // If we specify a fixed dt which is still >> fps, we can get
    // really good accuracies on interpolations of the midpoint
    //auto res2 = pr.doStep(res1.s, res1.et, res1.dt_next);
    // Integation dt of 10s still gives an accuracy of interpolation
    // of < 1.0E-5!!
    auto res2 = pr.doStep(res1.s, res1.et, 10);


    // do some interpolation in this interval and compare
    // with actual integrated values:
    PosState state_interp, state_integr;
    for(double f = 0.02; f < 1.0; f+= 0.02)
    {
        TimeDelta dtf = TimeDelta((res2.et - res1.et).value * f);
        //std::cout << dtf.value << std::endl;
        EphemerisTime etf = res1.et + dtf;

        // The interpolated state
        hermite(res1.s, res1.et, res2.s, res2.et, etf, state_interp);
        //std::cout << state_interp << std::endl;

        // the integrated state
        auto resf = pr.doStep(res1.s, res1.et, dtf);
        state_integr = resf.s;
        //std::cout << resf.s << std::endl;

        ASSERT_EQ(etf, resf.et);
        ASSERT_LT((state_integr.r - state_interp.r).length(), 1.0E-5);
        ASSERT_LT((state_integr.v - state_interp.v).length(), 1.0E-5);
    }

  
}




