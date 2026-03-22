
#include "../astro/ODE.h"
#include "../astro/State.h"
#include "../astro/Util.h"
#include <gtest/gtest.h>

#include <cmath>

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
    astro::Attractor earth = {Vec3(0.0), mu_earth};
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
    Mat3 I_singular = Mat3(0.0);
    RotODE rode;
    ASSERT_THROW(rode.setInertialMatrix(I_singular), astro::AstroException);

    Mat3 I_id = Mat3(1.0);
    ASSERT_NO_THROW(rode.setInertialMatrix(I_id));

    // Test inertia set in construction of rot ode
    ASSERT_THROW(RotODE rode2(I_singular), AstroException);
    ASSERT_NO_THROW(RotODE rode2(I_id));




}


