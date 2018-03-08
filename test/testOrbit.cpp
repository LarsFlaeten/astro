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


};



OrbitTest::OrbitTest()
  
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

TEST_F(OrbitTest, EmptyTest)
{
    ASSERT_EQ(1,1);
}




