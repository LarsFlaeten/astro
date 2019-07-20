#include "../astro/Util.h"
#include "../astro/State.cpp"
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
    s1.r = vec3d(1.0, 2.0, 3.0);
    s1.v = vec3d(4.0, 5.0, 6.0);
    PosState s2;
    s2.r = vec3d(1.0, 2.0, 3.0);
    s2.v = vec3d(4.0, 5.0, 6.0);
     
    PosState s3 = s1 + s2;

    ASSERT_EQ(s3.r, vec3d(2.0, 4.0, 6.0));
    ASSERT_EQ(s3.v, vec3d(8.0, 10.0, 12.0));
    
    PosState s4;
    s4.r = 0.5*s2.r;
    s4.v = 0.5*s2.v;

    PosState sn = s2 - s4;
    ASSERT_EQ(sn.r, vec3d(0.5, 1.0, 1.5));
    ASSERT_EQ(sn.v, vec3d(2.0, 2.5, 3.0));

    sn = 0.5*s2;
    sn = s2 - sn;
    ASSERT_EQ(sn.r, vec3d(0.5, 1.0, 1.5));
    ASSERT_EQ(sn.v, vec3d(2.0, 2.5, 3.0));

    sn = s2 + s4;
    ASSERT_EQ(sn.r, vec3d(1.5, 3.0, 4.5));
    ASSERT_EQ(sn.v, vec3d(6.0, 7.5, 9.0));


    // Divisor:
    sn = s4 / s1;
    ASSERT_EQ(sn.r, vec3d(0.5, 0.5, 0.5));
    ASSERT_EQ(sn.v, vec3d(0.5, 0.5, 0.5));


    // Abs
    sn = PosState(vec3d(-1, -2, -3), vec3d(-4, -5, -6));
    sn = abs(sn);
    ASSERT_EQ(sn.r, s1.r);
    ASSERT_EQ(sn.v, s1.v);

}

TEST_F(StateTest, TranslationInfNormTest)
{
    boost::numeric::odeint::vector_space_norm_inf< astro::PosState > norm_inf;
    PosState s1;
    s1.r = vec3d(1.0, 2.0, 3.0);
    s1.v = vec3d(4.0, 5.0, 6.0);

    
    ASSERT_EQ(norm_inf(s1), 6);


    s1.v.z = -21;
    ASSERT_EQ(norm_inf(s1), 21);
 
}

TEST_F(StateTest, QuatMulTest)
{
    quatd p(1, 2, 3, 4); p = p.normalize();
    double p0 = p.w; vec3d P(p.x, p.y, p.z);

    quatd q(3, 1, 4, 2); q = q.normalize();
    double q0 = q.w; vec3d Q(q.x, q.y, q.z);

    // Ork quat mul
    quatd r1 = p*q;
     
    // Ref quat mul
    double r20 = p0*q0 - P.dotproduct(Q);
    vec3d R2(p0*Q + q0*P + P.crossProduct(Q));
    quatd r2(R2.x, R2.y, R2.z, r20);

    ASSERT_LT(fabs(r1.x- r2.x), 1.0E-10);
    ASSERT_LT(fabs(r1.y- r2.y), 1.0E-10);
    ASSERT_LT(fabs(r1.z- r2.z), 1.0E-10);
    ASSERT_LT(fabs(r1.w- r2.w), 1.0E-10);





    //std::cout << "Ork: " << r1 << std::endl;
    //std::cout << "alt: " << r2 << std::endl;

}

TEST_F(StateTest, RotQuatTest1)
{
    vec3d v1(0.1, 0.5, 0.7);

    quatd q0(1, 0, 0, 1);
    q0 = q0.normalize();

    vec3d v2 = q0 * v1;
    quatd q2 = q0 * quatd(v1.x, v1.y, v1.z, 0) * q0.inverse();

    ASSERT_LT(fabs(v2.x- q2.x), 1.0E-10);
    ASSERT_LT(fabs(v2.y- q2.y), 1.0E-10);
    ASSERT_LT(fabs(v2.z- q2.z), 1.0E-10);


    //std::cout << v2 << std::endl;
    //std::cout << q2 << std::endl;


}

TEST_F(StateTest, OrkQuatMulVecBenchmark)
{
    vec3d v1(0.1, 0.5, 0.7);

    quatd q0(1, 0, 0, 1);
    q0 = q0.normalize();

    for(int i = 0; i < 10000000; i++)
        vec3d v2 = q0 * v1;
}

// 56% execution time compared to above... Allmost twice as fast!
TEST_F(StateTest, AltQuatMulVecBenchmark)
{
    vec3d v1(0.1, 0.5, 0.7);

    quatd q0(1, 0, 0, 1);
    q0 = q0.normalize();
    quatd q0_inv = q0.inverse();
    for(int i = 0; i < 10000000; i++)
        quatd q2 = q0 * quatd(v1.x, v1.y, v1.z, 0) * q0_inv;
}

TEST_F(StateTest, RotateTest1)
{

    // Reference vectors:
    vec3d ex(1, 0, 0);
    vec3d ey(0, 1, 0);
    vec3d ez(0, 0, 1);

    RotState rs;
    rs.q = quatd(ey, -astro::PIHALF); // A ship turned about y axis - 90 deg
    // At periapsis (rp, 0, 0), this shoul equate to a Nadir+ orientation
    // with the roof of the ship turned to the primary
    ASSERT_LT((rs.q*ex - ez).length(), 1.0E-6); // Local x is global z
    ASSERT_LT((rs.q*ey - ey).length(), 1.0E-6); // local y is global y
    ASSERT_LT((rs.q*ez - -ex).length(), 1.0E-6);// Local z is neg global x

    // Local angular velocity about z: (turing local left towards prograde)
    vec3d wb = astro::PIHALF*ez;

    // Global angular velocity should then be negative about global x:
    quatd w = rs.q * quatd(wb.x, wb.y, wb.z, 0.0) * rs.q.inverse();
    rs.w = vec3d(w.x, w.y, w.z);
    ASSERT_LT((rs.w.normalize() - -ex).length(), 1.0E-6);

    //std::cout << "wb: " << wb << std::endl;
    //std::cout << "w:  " << rs.w << std::endl;

    // A ship @PE turned prograde:
    rs.q = quatd(ez, astro::PIHALF);

    // We want our roof towards the primary
    // This would then be a right roll about local x (negative):
    wb = -astro::PIHALF*ex;
    w = rs.q * quatd(wb.x, wb.y, wb.z, 0.0) * rs.q.inverse();
    rs.w = vec3d(w.x, w.y, w.z);
    // And shold give us a negative turn bout global y:
    ASSERT_LT((rs.w.normalize() - -ey).length(), 1.0E-6);

   
    
    
}
void assert_almost_eq(mork::vec3d v1, mork::vec3d v2, double tol);

TEST_F(StateTest, TransformTest1)
{

    ASSERT_NO_THROW(astro::Spice().loadKernel("../data/spice/lsk/naif0012.tls"));
    ASSERT_NO_THROW(astro::Spice().loadKernel("../data/spice/spk/de430.bsp"));
    ASSERT_NO_THROW(astro::Spice().loadKernel("../data/spice/pck/pck00010.tpc"));
    astro::EphemerisTime et = astro::EphemerisTime::fromString("2018-06-12 21:00 UTC");;

    astro::PosState stateShip;
    stateShip.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
    stateShip.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]

    astro::ReferenceFrame rfj2000 = astro::ReferenceFrame::createJ2000();
    astro::ReferenceFrame rfearth = astro::ReferenceFrame::createBodyFixedSpice(399);


    astro::PosState newState = stateShip.transform(rfj2000, rfearth, et);

    std::cout << "Old state (J2000): " << stateShip.r << std::endl;
    std::cout << "New state (earth): " << newState.r << std::endl;

    ASSERT_LT(fabs(stateShip.r.length()- newState.r.length()), 1.0E-10);
    
    astro::PosState newState2 = newState.transform(rfearth, rfj2000, et);
    std::cout << "Old state (j2000): " << newState2.r << std::endl;



}



   
    
    






