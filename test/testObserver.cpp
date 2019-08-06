#include "../astro/Observer.cpp"
#include <gtest/gtest.h>

#include <cmath>


using namespace astro;

class ObserverTest : public ::testing::Test {

protected:
    ObserverTest();

    virtual ~ObserverTest();

    // Code here will be called immediately after the constructor (right
    // before each test).
    virtual void SetUp();

    // Code here will be called immediately after each test (right
    // before the destructor).
    virtual void TearDown();


};



ObserverTest::ObserverTest()
{
       
}

ObserverTest::~ObserverTest()
{

}

void ObserverTest::SetUp()
{

}

void ObserverTest::TearDown()
{
}

TEST_F(ObserverTest, InitTest)
{
    astro::State stateShip;
    stateShip.P.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
    stateShip.P.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]
    stateShip.R = astro::RotState();
    astro::Observer(399, astro::ReferenceFrame::createJ2000(), stateShip);
    
}

TEST_F(ObserverTest, ChangeFrameTest)
{
    ASSERT_NO_THROW(astro::Spice().loadKernel("../data/spice/lsk/naif0012.tls"));
    ASSERT_NO_THROW(astro::Spice().loadKernel("../data/spice/spk/de430.bsp"));
    ASSERT_NO_THROW(astro::Spice().loadKernel("../data/spice/pck/pck00010.tpc"));
    
    astro::State stateShip;
    stateShip.P.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
    stateShip.P.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]
    stateShip.R = astro::RotState();

    astro::Observer obs(399, astro::ReferenceFrame::createJ2000(), stateShip);

    std::cout << "J2000:" << std::endl;   
    std::cout << obs.getState().P << std::endl; 
    std::cout << "v: " << obs.getState().P.v.length() << "km/s" << std::endl;
    
    obs.setReferenceFrame(astro::ReferenceFrame::createBodyFixedSpice(399), true, 0);

    std::cout << "IAU_EARTH:" << std::endl;   
    std::cout << obs.getState().P << std::endl;
    std::cout << "v: " << obs.getState().P.v.length() << "km/s" << std::endl;
    
    obs.setReferenceFrame(astro::ReferenceFrame::createJ2000(), true, 0);

    std::cout << "J2000:" << std::endl;   
    std::cout << obs.getState().P << std::endl;
    std::cout << "v: " << obs.getState().P.v.length() << "km/s" << std::endl;
 
    obs.setReferenceFrame(astro::ReferenceFrame::createBodyFixedSpice(3), true, 0);

    std::cout << "IAU_EARTH_BC:" << std::endl;   
    std::cout << obs.getState().P << std::endl;
    std::cout << "v: " << obs.getState().P.v.length() << "km/s" << std::endl;
 
    obs.setReferenceFrame(astro::ReferenceFrame::createBodyFixedSpice(399), true, 0);

    std::cout << "IAU_EARTH:" << std::endl;   
    std::cout << obs.getState().P << std::endl;
    std::cout << "v: " << obs.getState().P.v.length() << "km/s" << std::endl;
 
    obs.setReferenceFrame(astro::ReferenceFrame::createBodyFixedSpice(3), true, 0);

    std::cout << "IAU_EARTH_BC:" << std::endl;   
    std::cout << obs.getState().P << std::endl;
    std::cout << "v: " << obs.getState().P.v.length() << "km/s" << std::endl;
 
    obs.setReferenceFrame(astro::ReferenceFrame::createJ2000(), true, 0);

    std::cout << "J2000:" << std::endl;   
    std::cout << obs.getState().P << std::endl;
    std::cout << "v: " << obs.getState().P.v.length() << "km/s" << std::endl;
 
}


