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

void assert_almost_neq(mork::vec3d v1, mork::vec3d v2, double tol)
{

    ASSERT_GT(fabs(v1.x-v2.x), tol);
    ASSERT_GT(fabs(v1.y-v2.y), tol);
    ASSERT_GT(fabs(v1.z-v2.z), tol);

    return;
}

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







    // Check outside time bonds for de430
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

TEST_F(SpiceCoreTest, getObserverRelativeStateAndPosition1)
{
    ASSERT_NO_THROW(astro::Spice().loadKernel("../data/spice/lsk/naif0012.tls"));
    ASSERT_NO_THROW(astro::Spice().loadKernel("../data/spice/spk/de430.bsp"));

    astro::EphemerisTime et = astro::EphemerisTime::fromString("2018-06-12 23:00 UTC");;
    astro::PosState stateErMBC;

    // Position of Earth rel to Mars BC
    astro::Spice().getRelativeGeometricState(399, 4, et, stateErMBC);

    astro::State stateShip;
    stateShip.P.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
    stateShip.P.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]

    astro::ReferenceFrame fr = astro::ReferenceFrame::createBodyFixedSpice(399);
     
    astro::ReferenceFrame fr2 = astro::ReferenceFrame::createBodyFixedSpice(3);
 
    astro::ReferenceFrame fr3 = astro::ReferenceFrame::createBodyFixedSpice(0);
   
    astro::ReferenceFrame fr4 = astro::ReferenceFrame::createJ2000();

    ASSERT_TRUE(fr4 == fr3);

    // Create observer @ earth in J2000 rotation frame
    astro::Observer obs(399, fr4, stateShip);    

    // Get state of Earth relative to Observer, should be inverse of stateShip:
    // No aberration correction
    astro::PosState earthRelState;
    astro::Spice().getRelativeState(399, obs, et, earthRelState);
    std::cout << "Earth state: " << earthRelState.r << std::endl;

    assert_almost_eq(stateShip.P.r, -earthRelState.r, 1.0E-10);
    assert_almost_eq(stateShip.P.v, -earthRelState.v, 1.0E-10);


    // Position of observer relative to Mars in J2000:
    astro::State stateShip2;
    stateShip2.P = earthRelState + stateErMBC;
    obs.setState(stateShip2);
    obs.setCenterObject(4);
    // Get state of Mars relative to Observer, should be inverse of stateShip:
    // No aberration correction
    astro::PosState marsRelState;
    astro::Spice().getRelativeState(4, obs, et, marsRelState);
    std::cout << "Mars state : " << marsRelState.r << std::endl;


    assert_almost_eq(stateShip2.P.r, -marsRelState.r, 1.0E-10);
    assert_almost_eq(stateShip2.P.v, -marsRelState.v, 1.0E-10);


}

// Check that a ship with the same position in different frames returns tha same state
TEST_F(SpiceCoreTest, getObserverRelativeStateAndPosition2)
{
    ASSERT_NO_THROW(astro::Spice().loadKernel("../data/spice/lsk/naif0012.tls"));
    ASSERT_NO_THROW(astro::Spice().loadKernel("../data/spice/spk/de430.bsp"));

    astro::EphemerisTime et = astro::EphemerisTime::fromString("2018-06-12 23:00 UTC");;
    astro::PosState stateErMBC;

    // Position of Earth rel to Mars BC
    astro::Spice().getRelativeGeometricState(399, 4, et, stateErMBC);

    astro::State stateShip;
    stateShip.P.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
    stateShip.P.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]


    astro::ReferenceFrame fr = astro::ReferenceFrame::createBodyFixedSpice(399);
     
    astro::ReferenceFrame fr2 = astro::ReferenceFrame::createBodyFixedSpice(3);
 
    astro::ReferenceFrame fr3 = astro::ReferenceFrame::createBodyFixedSpice(0);
   
    astro::ReferenceFrame fr4 = astro::ReferenceFrame::createJ2000();

    ASSERT_TRUE(fr4 == fr3);

    // Create observer @ earth in J2000 rotation frame
    astro::Observer obs(399, fr4, stateShip);    


    std::cout << "** No correction (geometric state) ** " << std::endl;
    // Get state of Mars relative to Observer:
    // No aberration correction
    astro::PosState marsRelState1;
    astro::Spice().getRelativeState(4, obs, et, marsRelState1);
    std::cout << "Mars state 1: " << marsRelState1.r << std::endl;


    // Position of observer relative to Mars in J2000:
    astro::State stateShip2;
    stateShip2.P = stateShip.P + stateErMBC;
    obs.setState(stateShip2);
    obs.setCenterObject(4);
    // Get state of Mars relative to Observer:
    // No aberration correction
    astro::PosState marsRelState2;
    astro::Spice().getRelativeState(4, obs, et, marsRelState2);
    std::cout << "Mars state 2: " << marsRelState2.r << std::endl;


    assert_almost_eq(marsRelState1.r, marsRelState2.r, 1.0E-10);
    assert_almost_eq(marsRelState1.v, marsRelState2.v, 1.0E-10);

    // Apply LT, the two states shoul still be equal:
    std::cout << "** LT correction ** " << std::endl;
 
    astro::Spice().getRelativeState(4, obs, et, marsRelState1, astro::AberrationCorrection::LightTime);
    std::cout << "Mars state 1: " << marsRelState1.r << std::endl;
    astro::Spice().getRelativeState(4, obs, et, marsRelState2, astro::AberrationCorrection::LightTime);
    std::cout << "Mars state 2: " << marsRelState2.r << std::endl;

    // Apply LT+S, the two states shoul still be equal:
    std::cout << "** LT+S correction (apparent state) ** " << std::endl;
 
    astro::Spice().getRelativeState(4, obs, et, marsRelState1, astro::AberrationCorrection::LightTimeStellar);
    std::cout << "Mars state 1: " << marsRelState1.r << std::endl;
    astro::Spice().getRelativeState(4, obs, et, marsRelState2, astro::AberrationCorrection::LightTimeStellar);
    std::cout << "Mars state 2: " << marsRelState2.r << std::endl;




}


// Check that a ship with the same position in different frames returns tha same state
TEST_F(SpiceCoreTest, getObserverRelativeStateAndPosition3)
{
    ASSERT_NO_THROW(astro::Spice().loadKernel("../data/spice/lsk/naif0012.tls"));
    ASSERT_NO_THROW(astro::Spice().loadKernel("../data/spice/spk/de430.bsp"));
    ASSERT_NO_THROW(astro::Spice().loadKernel("../data/spice/pck/pck00010.tpc"));
    astro::EphemerisTime et = astro::EphemerisTime::fromString("2018-06-12 23:00 UTC");;


    astro::State stateShip;
    stateShip.P.r = vec3d(0, 0, 0);  //[km]
    stateShip.P.v = vec3d(0, 0, 0);      //[km/s]


    astro::ReferenceFrame fr = astro::ReferenceFrame::createBodyFixedSpice(399);
     


    // Create observer @ earth in IAU_EARTH rotation frame
    astro::Observer obs(399, fr, stateShip);    


    // Get state of Mars relative to Observer (which is now earth):
    astro::PosState marsRelState1;
    astro::Spice().getRelativeState(4, obs, et, marsRelState1);

    // Use traditional function to retrive state between objects:
    mork::vec3d pos;
    astro::Spice().getRelativePosition(4, obs, et, pos);

    assert_almost_eq(pos, marsRelState1.r, 1.0E-10);     


    // TEst2:
    // Get state of Mars relative to Observer in J2000:
    // Check that state is different, but same length
    astro::PosState marsRelState2;
    obs.setReferenceFrame(astro::ReferenceFrame::createJ2000());
    astro::Spice().getRelativeState(4, obs, et, marsRelState2);
    
    assert_almost_neq(marsRelState1.r, marsRelState2.r, 1.0E-10);
    ASSERT_LT(fabs(marsRelState1.r.length() - marsRelState2.r.length()), 1.0E-10); 

    std::cout << "Mars state (J2000)     : " << marsRelState2.r << std::endl;
    std::cout << "Mars state (IAU_EARTH) : " << marsRelState1.r << std::endl;
 
}


TEST_F(SpiceCoreTest, getObserverRelativeStateCorrections)
{
    ASSERT_NO_THROW(astro::Spice().loadKernel("../data/spice/lsk/naif0012.tls"));
    ASSERT_NO_THROW(astro::Spice().loadKernel("../data/spice/spk/de430.bsp"));

    astro::EphemerisTime et = astro::EphemerisTime::fromString("2018-06-12 23:00 UTC");;

    astro::State stateShip;
    stateShip.P.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
    stateShip.P.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]


    astro::ReferenceFrame fr = astro::ReferenceFrame::createJ2000();

    // Create observer @ earth in J2000 rotation frame
    astro::Observer obs(399, fr, stateShip);    

    // Get state of Mars relative to Observer
    // No aberration correction
    astro::PosState marsStateNoC;
    astro::Spice().getRelativeState(4, obs, et, marsStateNoC);


    std::cout << "Mars, no correction  : " << marsStateNoC.r << std::endl;

    astro::PosState marsStateLT;
    astro::Spice().getRelativeState(4, obs, et, marsStateLT, astro::AberrationCorrection::LightTime );

    double A = marsStateNoC.r.length();
    double B = marsStateLT.r.length();
    double dotn = marsStateNoC.r.dotproduct(marsStateLT.r) / (A * B);
    std::cout << "Mars, LT correction  : " << marsStateLT.r << ", angle: " << acos(dotn)*1000 << " mRad (" << acos(dotn)*180.0*60*60/M_2_PI << " arc s)" << std::endl;

   astro::PosState marsStateLTS;
    astro::Spice().getRelativeState(4, obs, et, marsStateLTS, astro::AberrationCorrection::LightTimeStellar );

    double C = marsStateLTS.r.length();
    dotn = marsStateNoC.r.dotproduct(marsStateLTS.r) / (A * C);
 
    std::cout << "Mars, LTS correction : " << marsStateLTS.r << ", angle: " << acos(dotn)*1000 << " mRad (" << acos(dotn)*180.0*60*60/M_2_PI << " arc s)" << std::endl;



}

TEST_F(SpiceCoreTest, getPlanetaryConstantsTest)
{
    astro::Spice().loadKernel("../data/spice/pck/gm_de431.tpc");

    double gm;
    astro::Spice().getPlanetaryConstants(399, "GM", 1, &gm);
    ASSERT_LT(fabs(gm- 3.9860043543609598E+05), 1.0E-4);
    
    astro::Spice().getPlanetaryConstants(599, "GM", 1, &gm);
    ASSERT_LT(fabs(gm- 1.266865349218008E+08), 1.0E-4);
    
    double gm2;
    astro::Spice().getPlanetaryConstants(599, "GM", gm2);
    ASSERT_EQ(gm, gm2);
 
    astro::Spice().loadKernel("../data/spice/pck/pck00010.tpc");
    double radii[3];
    astro::Spice().getPlanetaryConstants(399, "RADII", 3, radii);
    ASSERT_GT(radii[0], 0.0);   
    ASSERT_GT(radii[1], 0.0);   
    ASSERT_GT(radii[2], 0.0);   


    mork::vec3d radii2;
    astro::Spice().getPlanetaryConstants(399, "RADII", radii2);
    ASSERT_EQ(radii[0], radii2.x);   
    ASSERT_EQ(radii[1], radii2.y);   
    ASSERT_EQ(radii[2], radii2.z);   


 

}


