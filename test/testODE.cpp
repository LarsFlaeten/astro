#include "../astro/ODE.cpp"

#include <gtest/gtest.h>

#include <cmath>

using namespace boost::numeric::odeint;

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
  :  mu_earth(398600.0), ode0(mu_earth)
{
		
   	


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

TEST_F(ODETest, EmptyTest)
{


}

