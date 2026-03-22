#include "../astro/Util.h"
#include "../astro/SpiceCore.h"
#include "../astro/ReferenceFrame.h"
#include "../astro/State.h"
#include "../astro/Time.h"
#include <gtest/gtest.h>

#include <cmath>


using namespace astro;

class StateTest : public ::testing::Test {

protected:
    StateTest();

    virtual ~StateTest();

    // Code here will be called immediately after the constructor (right
    // before each test).
    virtual void SetUp();

    // Code here will be called immediately after each test (right
    // before the destructor).
    virtual void TearDown();
};



StateTest::StateTest()
{

}

StateTest::~StateTest()
{

}

void StateTest::SetUp()
{
}

void StateTest::TearDown()
{
}

TEST_F(StateTest, TranslationOperatorTest1)
{
    PosState s1;
    s1.r = Vec3(1.0, 2.0, 3.0);
    s1.v = Vec3(4.0, 5.0, 6.0);
    PosState s2;
    s2.r = Vec3(1.0, 2.0, 3.0);
    s2.v = Vec3(4.0, 5.0, 6.0);
     
    PosState s3 = s1 + s2;

    ASSERT_EQ(s3.r, Vec3(2.0, 4.0, 6.0));
    ASSERT_EQ(s3.v, Vec3(8.0, 10.0, 12.0));
    
    PosState s4;
    s4.r = 0.5*s2.r;
    s4.v = 0.5*s2.v;

    PosState sn = s2 - s4;
    ASSERT_EQ(sn.r, Vec3(0.5, 1.0, 1.5));
    ASSERT_EQ(sn.v, Vec3(2.0, 2.5, 3.0));

    sn = 0.5*s2;
    sn = s2 - sn;
    ASSERT_EQ(sn.r, Vec3(0.5, 1.0, 1.5));
    ASSERT_EQ(sn.v, Vec3(2.0, 2.5, 3.0));

    sn = s2 + s4;
    ASSERT_EQ(sn.r, Vec3(1.5, 3.0, 4.5));
    ASSERT_EQ(sn.v, Vec3(6.0, 7.5, 9.0));


    // Divisor:
    sn = s4 / s1;
    ASSERT_EQ(sn.r, Vec3(0.5, 0.5, 0.5));
    ASSERT_EQ(sn.v, Vec3(0.5, 0.5, 0.5));


    // Abs
    sn = PosState(Vec3(-1, -2, -3), Vec3(-4, -5, -6));
    sn = abs(sn);
    ASSERT_EQ(sn.r, s1.r);
    ASSERT_EQ(sn.v, s1.v);

}

TEST_F(StateTest, QuatMulTest)
{
    Quat p(4, 1, 2, 3); p = glm::normalize(p);
    double p0 = p.w; Vec3 P(p.x, p.y, p.z);

    Quat q(2, 3, 1, 4); q = glm::normalize(q);
    double q0 = q.w; Vec3 Q(q.x, q.y, q.z);

    // GLM quat mul
    Quat r1 = p*q;

    // Ref quat mul
    double r20 = p0*q0 - glm::dot(P, Q);
    Vec3 R2(p0*Q + q0*P + glm::cross(P, Q));
    Quat r2(r20, R2.x, R2.y, R2.z);

    ASSERT_LT(fabs(r1.x- r2.x), 1.0E-10);
    ASSERT_LT(fabs(r1.y- r2.y), 1.0E-10);
    ASSERT_LT(fabs(r1.z- r2.z), 1.0E-10);
    ASSERT_LT(fabs(r1.w- r2.w), 1.0E-10);





    //std::cout << "Ork: " << r1 << std::endl;
    //std::cout << "alt: " << r2 << std::endl;

}

TEST_F(StateTest, RotQuatTest1)
{
    Vec3 v1(0.1, 0.5, 0.7);

    Quat q0(1, 1, 0, 0);
    q0 = glm::normalize(q0);

    Vec3 v2 = q0 * v1;
    Quat q2 = q0 * Quat(0, v1.x, v1.y, v1.z) * glm::inverse(q0);

    ASSERT_LT(fabs(v2.x- q2.x), 1.0E-10);
    ASSERT_LT(fabs(v2.y- q2.y), 1.0E-10);
    ASSERT_LT(fabs(v2.z- q2.z), 1.0E-10);


    //std::cout << v2 << std::endl;
    //std::cout << q2 << std::endl;


}

TEST_F(StateTest, RotateTest1)
{

    // Reference vectors:
    Vec3 ex(1, 0, 0);
    Vec3 ey(0, 1, 0);
    Vec3 ez(0, 0, 1);

    RotState rs;
    rs.q = glm::angleAxis(-astro::PIHALF, ey); // A ship turned about y axis - 90 deg
    // At periapsis (rp, 0, 0), this shoul equate to a Nadir+ orientation
    // with the roof of the ship turned to the primary
    ASSERT_LT(glm::length(rs.q*ex - ez), 1.0E-6); // Local x is global z
    ASSERT_LT(glm::length(rs.q*ey - ey), 1.0E-6); // local y is global y
    ASSERT_LT(glm::length(rs.q*ez - -ex), 1.0E-6);// Local z is neg global x

    // Local angular velocity about z: (turing local left towards prograde)
    Vec3 wb = astro::PIHALF*ez;

    // Global angular velocity should then be negative about global x:
    Quat w = rs.q * Quat(0.0, wb.x, wb.y, wb.z) * glm::inverse(rs.q);
    rs.w = Vec3(w.x, w.y, w.z);
    ASSERT_LT(glm::length(glm::normalize(rs.w) - -ex), 1.0E-6);

    //std::cout << "wb: " << wb << std::endl;
    //std::cout << "w:  " << rs.w << std::endl;

    // A ship @PE turned prograde:
    rs.q = glm::angleAxis(astro::PIHALF, ez);

    // We want our roof towards the primary
    // This would then be a right roll about local x (negative):
    wb = -astro::PIHALF*ex;
    w = rs.q * Quat(0.0, wb.x, wb.y, wb.z) * glm::inverse(rs.q);
    rs.w = Vec3(w.x, w.y, w.z);
    // And shold give us a negative turn bout global y:
    ASSERT_LT(glm::length(glm::normalize(rs.w) - -ey), 1.0E-6);

   
    
    
}
void assert_almost_eq(astro::Vec3 v1, astro::Vec3 v2, double tol);

TEST_F(StateTest, TransformTest1)
{

    ASSERT_NO_THROW(astro::Spice().loadKernel("../data/spice/lsk/naif0012.tls"));
    ASSERT_NO_THROW(astro::Spice().loadKernel("../data/spice/spk/de430.bsp"));
    ASSERT_NO_THROW(astro::Spice().loadKernel("../data/spice/pck/pck00010.tpc"));
    astro::EphemerisTime et = astro::EphemerisTime::fromString("2018-06-12 21:00 UTC");;

    astro::State stateShip;
    stateShip.P.r = Vec3(-6045.0, -3490.0, 2500.0);  //[km]
    stateShip.P.v = Vec3(-3.457, 6.618, 2.533);      //[km/s]

    astro::ReferenceFrame rfj2000 = astro::ReferenceFrame::createJ2000();
    astro::ReferenceFrame rfearth = astro::ReferenceFrame::createBodyFixedSpice(399);


    astro::State newState = stateShip.transform(rfj2000, rfearth, et);

    std::cout << "Old state (J2000): " << stateShip.P.r << std::endl;
    std::cout << "New state (earth): " << newState.P.r << std::endl;

    ASSERT_LT(fabs(glm::length(stateShip.P.r)- glm::length(newState.P.r)), 1.0E-10);
    
    astro::State newState2 = newState.transform(rfearth, rfj2000, et);
    std::cout << "Old state (j2000): " << newState2.P.r << std::endl;



}



   
    
    






