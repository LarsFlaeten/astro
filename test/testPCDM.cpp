#include "../astro/State.h"
#include "../astro/Time.h"
#include "../astro/PCDM.cpp"

#include <gtest/gtest.h>

#include <cmath>

using namespace boost::numeric::odeint;
using namespace astro;


// This test is only testing the interface of PCDM
// Actual numerical integrations are tested through
// propagator in NumIntTest
class PCDMTest : public ::testing::Test {

protected:
    PCDMTest();

    virtual ~PCDMTest();

    // Code here will be called immediately after the constructor (right
    // before each test).
    virtual void SetUp();

    // Code here will be called immediately after each test (right
    // before the destructor).
    virtual void TearDown();

};



PCDMTest::PCDMTest()
{
    


}

PCDMTest::~PCDMTest()
{

}

void PCDMTest::SetUp()
{
}

void PCDMTest::TearDown()
{
}



TEST_F(PCDMTest, PCDMTest)
{

    // No test so far    
}


