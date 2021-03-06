#include "../astro/State.h"
#include "../astro/Time.h"
#include "../astro/Propagator.h"
#include "../astro/RKF78.cpp"
#include "../astro/RKF45.cpp"
#include "../astro/Orbit.h"
#include "../astro/ODE.h"
#include "../astro/PCDM.h"

#include <gtest/gtest.h>

#include <cmath>

using namespace boost::numeric::odeint;
using namespace astro;

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

    astro::PosState    state0;
    astro::OrbitElements oe0;
    astro::EphemerisTime et0;
    double mu_earth;
    astro::Attractor earth;
    astro::ODE  ode0;
};



NumIntTest::NumIntTest()
  :  et0(0), mu_earth(398600.0)
{
    earth = {vec3d::ZERO, mu_earth};
    ode0.addAttractor(earth);

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
    
    astro::Propagator<astro::ODE, astro::RKF78> pr(ode0);
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
    
    astro::Propagator<astro::ODE, astro::RKF45> pr(ode0);
    astro::RKF45::setTolerance(1.0E-8);       

    auto resv = pr.doSteps(state0, et0, et1, astro::TimeDelta(DT));
}

TEST_F(NumIntTest, RKF78HypAsymptote)
{
    astro::PosState s;
    s.r = vec3d(6378 + 400, 0.0, 0.0);
    //double v = 7.66895; // Circular
    double v = 10.84509; // roughly eccentric (to six decimals on e)
    s.v = vec3d(0, 0, v*1.1);
    

    astro::OrbitElements oe = astro::OrbitElements::fromStateVectorOE(s, et0, mu_earth);
    //std::cout << oe << std::endl; 

    // Propagate the orbit for a long time (12 days), far outside SOI:
    astro::EphemerisTime et1 = et0 + 12*60*60*24;

    astro::Propagator<astro::ODE, astro::RKF78> pr(ode0);
    astro::RKF78::setTolerance(1.0E-8);       

    double DT = 1.0; // Initial dt
    auto resv = pr.doSteps(s, et0, et1, astro::TimeDelta(DT));
    
    //for(auto res : resv)
    //    std::cout << res.dt_next.value << "\t" << res.s.r.length() << "\t" << res.s.v.length() << std::endl;

    astro::PosState s_inf = resv.back().s;

    double v_inf = s_inf.v.length();
    ASSERT_LT(fabs(v_inf- astro::hyperbolicExcessVelocity(mu_earth, oe.a)), 0.016);
    
    // asymptotic angles are more easy based on velocity (fast convergence)
    // Position based asymptote is much more slow-convergent (about same as v_inf convergence above)
    double hyp_asym = astro::hyperbolicAsymptote(oe.e);
    vec3d v0 = -s.r.normalize();
    vec3d v1 = s_inf.v.normalize();
    vec3d v2 = s_inf.r.normalize();
    double hyp_asym1 = acos(v0.dotproduct(v1)); // Asymptote from velocity
    double hyp_asym2 = acos(v0.dotproduct(v2)); // Asymptote from position
  
    ASSERT_LT(fabs(hyp_asym - hyp_asym1), 0.0005); // 0.03 degrees..
    ASSERT_LT(fabs(hyp_asym - hyp_asym2), 0.0032); // 0.2 degrees..
}


void assert_almost_eq(mork::vec3d v1, mork::vec3d v2, double tol)
{

    ASSERT_LT(fabs(v1.x-v2.x), tol);
    ASSERT_LT(fabs(v1.y-v2.y), tol);
    ASSERT_LT(fabs(v1.z-v2.z), tol);

    return;
}
// Constant angular velocity test
TEST_F(NumIntTest, PCDMTestConstW)
{
        double w = 0.1;
        astro::RotState rs;
        rs.q = quatd(0, 0, 0, 1); // "Unit" quaternion

        // The differential equation for rotations:
        astro::RotODE rode;

        astro::EphemerisTime et(12345);
        // Find the time when we shold be rotated 90 degrees by the axis:
        double T = astro::PIHALF / w;
        astro::EphemerisTime et2 = et + TimeDelta(T);
        

        astro::TimeDelta dt(1/60.0);
        //astro::TimeDelta dt(1.0/60);

        // ---
        // Rotation about positive global X
        // ---
        rs.w = vec3d(w, 0.0, 0.0); // rotation about global X by 0.1 rad/s
        
        // Do the integration 
        auto resv = astro::PCDM::doSteps(rode, rs, et, et2, dt);

        // After this, x will be x, y will be z, and z will be -y.
        // Body vectors
        vec3d ex = vec3d::UNIT_X;
        vec3d ey = vec3d::UNIT_Y;
        vec3d ez = vec3d::UNIT_Z;

        // Transformt to global
        ex = resv.back().rs.q * ex;
        ey = resv.back().rs.q * ey;
        ez = resv.back().rs.q * ez;

        assert_almost_eq(ex, vec3d::UNIT_X, 1.0E-10);
        assert_almost_eq(ey, vec3d::UNIT_Z, 1.0E-10);
        assert_almost_eq(ez, -1.0*vec3d::UNIT_Y, 1.0E-10);

        // Lenght of quat should be 1
        ASSERT_LT(fabs(resv.back().rs.q.length()-1.0), 1.0E-10);

        // ---
        // Rotation about negative global X
        // ---
        rs.w = vec3d(-w, 0.0, 0.0); // rotation about global X by 0.1 rad/s
        
        // Do the integration 
        resv = astro::PCDM::doSteps(rode, rs, et, et2, dt);

        // After this, x will be x, y will be -z, and z will be y.
        // body vectors
        ex = vec3d::UNIT_X;
        ey = vec3d::UNIT_Y;
        ez = vec3d::UNIT_Z;

        // Transform to global:
        ex = resv.back().rs.q * ex;
        ey = resv.back().rs.q * ey;
        ez = resv.back().rs.q * ez;

        assert_almost_eq(ex, vec3d::UNIT_X, 1.0E-10);
        assert_almost_eq(ey, -1.0*vec3d::UNIT_Z, 1.0E-10);
        assert_almost_eq(ez, vec3d::UNIT_Y, 1.0E-10);

        // Lenght of quat should be 1
        ASSERT_LT(fabs(resv.back().rs.q.length()-1.0), 1.0E-10);

        // ---
        // Rotation about positive global y
        // ---
        rs.w = vec3d(0.0, w, 0.0);
        
        // Do the integration 
        resv = astro::PCDM::doSteps(rode, rs, et, et2, dt);

        // After this, x will be -z, y will be y, and z will be x.
        // body vectors
        ex = vec3d::UNIT_X;
        ey = vec3d::UNIT_Y;
        ez = vec3d::UNIT_Z;

        // Transform to global:
        ex = resv.back().rs.q * ex;
        ey = resv.back().rs.q * ey;
        ez = resv.back().rs.q * ez;

        assert_almost_eq(ex, -1.0*vec3d::UNIT_Z, 1.0E-10);
        assert_almost_eq(ey, vec3d::UNIT_Y, 1.0E-10);
        assert_almost_eq(ez, vec3d::UNIT_X, 1.0E-10);

        // Lenght of quat should be 1
        ASSERT_LT(fabs(resv.back().rs.q.length()-1.0), 1.0E-10);

        // ---
        // Rotation about negative global y
        // ---
        rs.w = vec3d(0.0, -w, 0.0);
        
        // Do the integration 
        resv = astro::PCDM::doSteps(rode, rs, et, et2, dt);

        // After this, x will be z, y will be y, and z will be -x.
        // body vectors
        ex = vec3d::UNIT_X;
        ey = vec3d::UNIT_Y;
        ez = vec3d::UNIT_Z;

        // Transform to global:
        ex = resv.back().rs.q * ex;
        ey = resv.back().rs.q * ey;
        ez = resv.back().rs.q * ez;

        assert_almost_eq(ex, vec3d::UNIT_Z, 1.0E-10);
        assert_almost_eq(ey, vec3d::UNIT_Y, 1.0E-10);
        assert_almost_eq(ez, -1.0*vec3d::UNIT_X, 1.0E-10);

        // Lenght of quat should be 1
        ASSERT_LT(fabs(resv.back().rs.q.length()-1.0), 1.0E-10);

        // ---
        // Rotation about positive global z
        // ---
        rs.w = vec3d(0.0, 0.0, w);
        
        // Do the integration 
        resv = astro::PCDM::doSteps(rode, rs, et, et2, dt);

        // After this, x will be y, y will be -x, and z will be z.
        // body vectors
        ex = vec3d::UNIT_X;
        ey = vec3d::UNIT_Y;
        ez = vec3d::UNIT_Z;

        // Transform to global:
        ex = resv.back().rs.q * ex;
        ey = resv.back().rs.q * ey;
        ez = resv.back().rs.q * ez;

        assert_almost_eq(ex, vec3d::UNIT_Y, 1.0E-10);
        assert_almost_eq(ey, -1.0*vec3d::UNIT_X, 1.0E-10);
        assert_almost_eq(ez, vec3d::UNIT_Z, 1.0E-10);

        // Lenght of quat should be 1
        ASSERT_LT(fabs(resv.back().rs.q.length()-1.0), 1.0E-10);

        // ---
        // Rotation about negative global z
        // ---
        rs.w = vec3d(0.0, 0.0, -w);
        
        // Do the integration 
        resv = astro::PCDM::doSteps(rode, rs, et, et2, dt);

        // After this, x will be -y, y will be x, and z will be z.
        // body vectors
        ex = vec3d::UNIT_X;
        ey = vec3d::UNIT_Y;
        ez = vec3d::UNIT_Z;

        // Transform to global:
        ex = resv.back().rs.q * ex;
        ey = resv.back().rs.q * ey;
        ez = resv.back().rs.q * ez;

        assert_almost_eq(ex, -1.0*vec3d::UNIT_Y, 1.0E-10);
        assert_almost_eq(ey, vec3d::UNIT_X, 1.0E-10);
        assert_almost_eq(ez, vec3d::UNIT_Z, 1.0E-10);

        // Lenght of quat should be 1
        ASSERT_LT(fabs(resv.back().rs.q.length()-1.0), 1.0E-10);
}


// Constant angular velocity test
TEST_F(NumIntTest, PCDMTestZeroW)
{
        double w = 0.1;
        astro::RotState rs;
        rs.q = quatd(0, 0, 0, 1); // "Unit" quaternion

        // The differential equation for rotations:
        astro::RotODE rode;

        astro::EphemerisTime et(12345);
        // Find the time when we shold be rotated 90 degrees by the axis:
        double T = astro::PIHALF / w;
        astro::EphemerisTime et2 = et + TimeDelta(T);
        

        astro::TimeDelta dt(1/60.0);
        //astro::TimeDelta dt(1.0/60);


        // ---
        // No Rotation should retain orientation
        // ---
        rs.w = vec3d(0.0, 0.0, 0.0);
        
        // Do the integration 
        auto resv = astro::PCDM::doSteps(rode, rs, et, et2, dt);

        // After this, x will be x, y will be y, and z will be z.
        // body vectors
        vec3d ex = vec3d::UNIT_X;
        vec3d ey = vec3d::UNIT_Y;
        vec3d ez = vec3d::UNIT_Z;

        std::cout << resv.back().rs.q << std::endl;

        // Transform to global:
        ex = resv.back().rs.q * ex;
        ey = resv.back().rs.q * ey;
        ez = resv.back().rs.q * ez;

        assert_almost_eq(ex, vec3d::UNIT_X, 1.0E-10);
        assert_almost_eq(ey, vec3d::UNIT_Y, 1.0E-10);
        assert_almost_eq(ez, vec3d::UNIT_Z, 1.0E-10);

        // Lenght of quat should be 1
        ASSERT_LT(fabs(resv.back().rs.q.length()-1.0), 1.0E-10);


}

// Application of torque positive x test
TEST_F(NumIntTest, PCDMTestTorqueXP)
{
        double w = 0.1; // The average angular velocity we want
        astro::RotState rs;
        rs.q = quatd(0, 0, 0, 1); // "Unit" quaternion
        rs.w = vec3d(0.0, 0.0, 0.0); // zero angular velocity
 
        // The differential equation for rotations:
        astro::RotODE rode;

        astro::EphemerisTime et(12345);
        // Find the time when we shold be rotated 90 degrees by the axis:
        double T = astro::PIHALF / w;
        astro::EphemerisTime et2 = et + TimeDelta(T);
        
        //Apply torque
        rode.setBodyTorque(vec3d(2*w/T, 0, 0));

       

        astro::TimeDelta dt(1/60.0);

        // Do the integration 
        auto resv = astro::PCDM::doSteps(rode, rs, et, et2, dt);


        // After this, x will be x, y will be z, and z will be -y.
        // body vectors
        vec3d ex = vec3d::UNIT_X;
        vec3d ey = vec3d::UNIT_Y;
        vec3d ez = vec3d::UNIT_Z;

        // Transform to global:
        ex = resv.back().rs.q * ex;
        ey = resv.back().rs.q * ey;
        ez = resv.back().rs.q * ez;

        assert_almost_eq(ex, vec3d::UNIT_X, 1.0E-10);
        assert_almost_eq(ey, vec3d::UNIT_Z, 1.0E-10);
        assert_almost_eq(ez, -vec3d::UNIT_Y, 1.0E-10);
}
 
// Application of torque negative x test
TEST_F(NumIntTest, PCDMTestTorqueXN)
{
        double w = 0.1; // The average angular velocity we want
        astro::RotState rs;
        rs.q = quatd(0, 0, 0, 1); // "Unit" quaternion
        rs.w = vec3d(0.0, 0.0, 0.0); // zero angular velocity
 
        // The differential equation for rotations:
        astro::RotODE rode;

        astro::EphemerisTime et(12345);
        // Find the time when we shold be rotated 90 degrees by the axis:
        double T = astro::PIHALF / w;
        astro::EphemerisTime et2 = et + TimeDelta(T);
        
        //Apply torque
        rode.setBodyTorque(vec3d(-2*w/T, 0, 0));

       

        astro::TimeDelta dt(1/60.0);

        // Do the integration 
        auto resv = astro::PCDM::doSteps(rode, rs, et, et2, dt);


        // After this, x will be x, y will be -z, and z will be y.
        // body vectors
        vec3d ex = vec3d::UNIT_X;
        vec3d ey = vec3d::UNIT_Y;
        vec3d ez = vec3d::UNIT_Z;

        // Transform to global:
        ex = resv.back().rs.q * ex;
        ey = resv.back().rs.q * ey;
        ez = resv.back().rs.q * ez;

        assert_almost_eq(ex, vec3d::UNIT_X, 1.0E-10);
        assert_almost_eq(ey, -vec3d::UNIT_Z, 1.0E-10);
        assert_almost_eq(ez, vec3d::UNIT_Y, 1.0E-10);
}
 
