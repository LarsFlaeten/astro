#include "../astro/SpiceCore.cpp"
#include <gtest/gtest.h>

class SpiceCoreTest : public ::testing::Test {

protected:
    SpiceCoreTest();

    virtual ~SpiceCoreTest();

    // Code here will be called immediately after the constructor (right
    // before each test).
    virtual void SetUp();

    // Code here will be called immediately after each test (right
    // before the destructor).
    virtual void TearDown();
};



SpiceCoreTest::SpiceCoreTest()
{

}

SpiceCoreTest::~SpiceCoreTest()
{

}

void SpiceCoreTest::SetUp()
{

}

void SpiceCoreTest::TearDown()
{

}

TEST_F(SpiceCoreTest, LoadKernelTest1)
{

    ASSERT_NO_THROW(astro::Spice().loadKernel("../data/spice/lsk/naif0012.tls"));

    ASSERT_THROW(astro::Spice().loadKernel("NoFileWithThisName"), astro::SpiceException);    
}
