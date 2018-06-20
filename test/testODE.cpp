#include "../astro/ODE.cpp"

#include <gtest/gtest.h>

#include <cmath>

using namespace boost::numeric::odeint;
using namespace astro;

class ODETest : public ::testing::Test {

protected:
    ODETest();

    virtual ~ODETest();

    // Code here will be called immediately after the constructor (right
    // before each test).
    virtual void SetUp();

    // Code here will be called immediately after each test (right
    // before the destructor).
    virtual void TearDown();

    double mu_earth;
    astro::ODE  ode0;
};



ODETest::ODETest()
  :  mu_earth(398600.0)
{
    astro::Attractor earth = {vec3d::ZERO, mu_earth};
    ode0.addAttractor(earth);
   	


}

ODETest::~ODETest()
{

}

void ODETest::SetUp()
{
}

void ODETest::TearDown()
{
}

TEST_F(ODETest, RotODESingularInertia)
{
    ork::mat3d I_singular = ork::mat3d::ZERO;
    RotODE rode;
    ASSERT_THROW(rode.setInertialMatrix(I_singular), astro::AstroException);

    ork::mat3d I_id = ork::mat3d::IDENTITY;
    ASSERT_NO_THROW(rode.setInertialMatrix(I_id));

    // Test inertia set in construction of rot ode
    ASSERT_THROW(RotODE rode2(I_singular), AstroException);
    ASSERT_NO_THROW(RotODE rode2(I_id));



    
}

