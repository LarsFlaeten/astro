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

    // Test loading same kernel twice, Spice should handle this fine
    ASSERT_NO_THROW(astro::Spice().loadKernel("../data/spice/lsk/naif0012.tls"));
}

void assert_almost_eq(ork::vec3d v1, ork::vec3d v2, double tol);

TEST_F(SpiceCoreTest, getStateTest1)
{
    ASSERT_NO_THROW(astro::Spice().loadKernel("../data/spice/lsk/naif0012.tls"));
    ASSERT_NO_THROW(astro::Spice().loadKernel("../data/spice/spk/de430.bsp"));

    astro::EphemerisTime et = astro::EphemerisTime::fromString("2018-06-12 23:00 UTC");;
    std::cout << et.toISOUTCString() << std::endl;
    astro::PosState state1, state2;
    ASSERT_NO_THROW(astro::Spice().getRelativeGeometricState(399,0, et, state1));  

    std::cout << "position of Earth in solar system:\n" << state1 << std::endl;

    // Advance 1s:
    et += astro::TimeDelta(1.0);
    ASSERT_NO_THROW(astro::Spice().getRelativeGeometricState(399,0, et, state2));  

    ork::vec3d diff = state2.r - state1.r;
    assert_almost_eq(diff, state1.v, 1.0E-4);


    // Check outside tome bonds for de430
    et = astro::EphemerisTime::fromString("2650 Jan 25");
    ASSERT_THROW(astro::Spice().getRelativeGeometricState(399,0,et,state1), astro::SpiceException);


}

TEST_F(SpiceCoreTest, getGeometricStateBenchmarkTestEarthSSB1)
{
    // Retreive 100000 states:
    astro::EphemerisTime et = astro::EphemerisTime::fromString("2018-06-12 23:00 UTC");;
    astro::PosState state1, state2;

    for(int i = 0; i < 100000; ++i)
        astro::Spice().getRelativeGeometricState(399,0, et, state1);  

   



}
TEST_F(SpiceCoreTest, getGeometricStateBenchmarkTestEarthMoon)
{
    // Retreive 100000 states:
    astro::EphemerisTime et = astro::EphemerisTime::fromString("2018-06-12 23:00 UTC");;
    astro::PosState state1, state2;

    for(int i = 0; i < 100000; ++i)
        astro::Spice().getRelativeGeometricState(399,301, et, state1);  

   



}


TEST_F(SpiceCoreTest, getGeometricStateBenchmarkTestEarthSSB2)
{
    // Retreive 100000 states:
    astro::EphemerisTime et = astro::EphemerisTime::fromString("2018-06-12 23:00 UTC");;
    astro::PosState state1, state2;

    for(int i = 0; i < 100000; ++i) {
        et += astro::TimeDelta(1.0);
        astro::Spice().getRelativeGeometricState(399,0, et, state1);  
    }
   



}
