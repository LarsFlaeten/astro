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
