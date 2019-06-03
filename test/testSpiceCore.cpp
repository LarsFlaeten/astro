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

void assert_almost_eq(mork::vec3d v1, mork::vec3d v2, double tol);

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

    mork::vec3d diff = state2.r - state1.r;
    assert_almost_eq(diff, state1.v, 1.0E-4);


    // Check outside tome bonds for de430
    et = astro::EphemerisTime::fromString("2650 Jan 25");
    ASSERT_THROW(astro::Spice().getRelativeGeometricState(399,0,et,state1), astro::SpiceException);


}

TEST_F(SpiceCoreTest, debugDumpTest)
{
    ASSERT_NO_THROW(astro::Spice().loadKernel("../data/spice/lsk/naif0012.tls"));
    ASSERT_NO_THROW(astro::Spice().loadKernel("../data/spice/spk/de430.bsp"));

    std::cout << astro::Spice() << std::endl;
    
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

TEST_F(SpiceCoreTest, getGeometricStateEarthMars)
{
    // Retreive 100000 states:
    astro::EphemerisTime et = astro::EphemerisTime::fromString("2018-06-12 23:00 UTC");;
    astro::PosState state1, state2, state3;
    
    // Position of earth reletive to SSB
    astro::Spice().getRelativeGeometricState(399,0, et, state1);  
    // Position of mars relative to ssb
    astro::Spice().getRelativeGeometricState(4,0, et, state2);

    // Position of earth as seen from mars
    astro::Spice().getRelativeGeometricState(399, 4, et, state3);
    
    assert_almost_eq(state1.r - state2.r, state3.r, 1.0E-6); 
    assert_almost_eq(state1.v - state2.v, state3.v, 1.0E-6);

}

TEST_F(SpiceCoreTest, getPositionTest1)
{
    ASSERT_NO_THROW(astro::Spice().loadKernel("../data/spice/lsk/naif0012.tls"));
    ASSERT_NO_THROW(astro::Spice().loadKernel("../data/spice/spk/de430.bsp"));

    astro::EphemerisTime et = astro::EphemerisTime::fromString("2018-06-12 23:00 UTC");;
    astro::PosState state1, state2;

    // Position of earth from SSB
    ASSERT_NO_THROW(astro::Spice().getRelativeGeometricState(399,0, et, state1));  
    vec3d pos1, pos2, pos3;
    ASSERT_NO_THROW(astro::Spice().getRelativePosition(399,0,et, pos1, astro::AberrationCorrection::None));
    assert_almost_eq(state1.r, pos1, 1.0E-10);

    // Position of earth from mars
    ASSERT_NO_THROW(astro::Spice().getRelativeGeometricState(399,4, et, state2));  
    ASSERT_NO_THROW(astro::Spice().getRelativePosition(399,4,et, pos1, astro::AberrationCorrection::None));
    ASSERT_NO_THROW(astro::Spice().getRelativePosition(399,4,et, pos2));
    ASSERT_NO_THROW(astro::Spice().getRelativePosition(399,4,et, pos3, astro::AberrationCorrection::LightTime));
    assert_almost_eq(state2.r, pos1, 1.0E-10);
    assert_almost_eq(pos1, pos2, 1.0E-10);


    // Pos3 should be offset from pos2, since it was calculated with light time correction
    ASSERT_NE(pos2.x, pos3.x);
    ASSERT_NE(pos2.y, pos3.y);
    ASSERT_NE(pos2.z, pos3.z);







    // Check outside tome bonds for de430
    et = astro::EphemerisTime::fromString("2650 Jan 25");
    ASSERT_THROW(astro::Spice().getRelativePosition(399,0,et,pos1), astro::SpiceException);


}

TEST_F(SpiceCoreTest, getPositionTestBenchMarkAbCorrNone)
{
    // Retreive 100000 states:
    astro::EphemerisTime et = astro::EphemerisTime::fromString("2018-06-12 23:00 UTC");;
    mork::vec3d pos;

    for(int i = 0; i < 100000; ++i) {
        et += astro::TimeDelta(1.0);
        astro::Spice().getRelativePosition(399,0, et, pos,
                astro::AberrationCorrection::None);  
    }
   
}

TEST_F(SpiceCoreTest, getPositionTestBenchMarkAbCorrLightTime)
{
    // Retreive 100000 states:
    astro::EphemerisTime et = astro::EphemerisTime::fromString("2018-06-12 23:00 UTC");;
    mork::vec3d pos;

    for(int i = 0; i < 100000; ++i) {
        et += astro::TimeDelta(1.0);
        astro::Spice().getRelativePosition(399,0, et, pos,
                astro::AberrationCorrection::LightTime);  
    }
   
}
TEST_F(SpiceCoreTest, getPositionTestBenchMarkAbCorrLightTimeStellar)
{
    // Retreive 100000 states:
    astro::EphemerisTime et = astro::EphemerisTime::fromString("2018-06-12 23:00 UTC");;
    mork::vec3d pos;

    for(int i = 0; i < 100000; ++i) {
        et += astro::TimeDelta(1.0);
        astro::Spice().getRelativePosition(399,0, et, pos,
                astro::AberrationCorrection::LightTimeStellar);  
    }
   
}

TEST_F(SpiceCoreTest, getPositionTestBenchMarkAbCorrCNLightTime)
{
    // Retreive 100000 states:
    astro::EphemerisTime et = astro::EphemerisTime::fromString("2018-06-12 23:00 UTC");;
    mork::vec3d pos;

    for(int i = 0; i < 100000; ++i) {
        et += astro::TimeDelta(1.0);
        astro::Spice().getRelativePosition(399,0, et, pos,
                astro::AberrationCorrection::CNLightTime);  
    }
   
}
TEST_F(SpiceCoreTest, getPositionTestBenchMarkAbCorrCNLightTimeStellar)
{
    // Retreive 100000 states:
    astro::EphemerisTime et = astro::EphemerisTime::fromString("2018-06-12 23:00 UTC");;
    mork::vec3d pos;

    for(int i = 0; i < 100000; ++i) {
        et += astro::TimeDelta(1.0);
        astro::Spice().getRelativePosition(399,0, et, pos,
                astro::AberrationCorrection::CNLightTimeStellar);  
    }
   
}

TEST_F(SpiceCoreTest, getPlanetaryConstantsTest)
{
    astro::Spice().loadKernel("../data/spice/pck/gm_de431.tpc");

    double gm;
    astro::Spice().getPlanetaryConstants(399, "GM", 1, &gm);
    ASSERT_LT(fabs(gm- 3.9860043543609598E+05), 1.0E-4);
    
    astro::Spice().getPlanetaryConstants(599, "GM", 1, &gm);
    ASSERT_LT(fabs(gm- 1.266865349218008E+08), 1.0E-4);
    
    astro::Spice().loadKernel("../data/spice/pck/pck00010.tpc");
    double radii[3];
    astro::Spice().getPlanetaryConstants(399, "RADII", 3, radii);
    ASSERT_GT(radii[0], 0.0);   
    ASSERT_GT(radii[1], 0.0);   
    ASSERT_GT(radii[2], 0.0);   
 
}


