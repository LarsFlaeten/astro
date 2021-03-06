#include "../astro/Util.h"
#include "../astro/ReferenceFrame.cpp"
#include <gtest/gtest.h>

#include <cmath>

#include <algorithm>
#include <string>

using namespace astro;

class ReferenceFrameTest : public ::testing::Test {

protected:
    ReferenceFrameTest();

    virtual ~ReferenceFrameTest();

    // Code here will be called immediately after the constructor (right
    // before each test).
    virtual void SetUp();

    // Code here will be called immediately after each test (right
    // before the destructor).
    virtual void TearDown();
};



ReferenceFrameTest::ReferenceFrameTest()
{

}

ReferenceFrameTest::~ReferenceFrameTest()
{

}

void ReferenceFrameTest::SetUp()
{
}

void ReferenceFrameTest::TearDown()
{
}


TEST_F(ReferenceFrameTest, InertialFrameTest1)
{
    // CHeck that two J2000 frames gives a relative orientaion that equals identity
    astro::ReferenceFrame rf = ReferenceFrame::createJ2000();
    astro::ReferenceFrame rf2 = ReferenceFrame::createJ2000();

    ASSERT_EQ(rf, rf2);
    mork::mat3d rrot = rf.getRotationToJ2000(EphemerisTime(0));
    ASSERT_EQ(rrot==mork::mat3d::IDENTITY, true);
    ASSERT_EQ(rf.getName(), "J2000");
    ASSERT_EQ(rf.getId(), 1);



}

TEST_F(ReferenceFrameTest, IsJ2000Test)
{
    astro::ReferenceFrame rf1 = ReferenceFrame::createJ2000();
    astro::ReferenceFrame rf2 = ReferenceFrame::createBodyFixedSpice(3);
    astro::ReferenceFrame rf3 = ReferenceFrame::createBodyFixedSpice(399);
    astro::ReferenceFrame rf4 = ReferenceFrame::createBodyFixedSpice(301);
    astro::ReferenceFrame rf5 = ReferenceFrame::createBodyFixedSpice(5);
    astro::ReferenceFrame rf6 = ReferenceFrame::createBodyFixedSpice(599);

    std::cout << "Frames:" << std::endl;
    std::cout << rf1.getName() << std::endl;
    std::cout << rf2.getName() << std::endl;
    std::cout << rf3.getName() << std::endl;
    std::cout << rf4.getName() << std::endl;
    std::cout << rf5.getName() << std::endl;
    std::cout << rf6.getName() << std::endl;

    ASSERT_TRUE(rf1.isJ2000());
    ASSERT_FALSE(rf2.isJ2000());
    ASSERT_FALSE(rf3.isJ2000());
    ASSERT_FALSE(rf4.isJ2000());
    ASSERT_FALSE(rf5.isJ2000());
    ASSERT_FALSE(rf6.isJ2000());
}


TEST_F(ReferenceFrameTest, InertialFrameTest2)
{
    // CHeck that a body fixed frame created arouund
    // body 0 (SSB) is the same as the J2000 frame.
    astro::ReferenceFrame rf = ReferenceFrame::createBodyFixedSpice(0);
    astro::ReferenceFrame rf2 = ReferenceFrame::createJ2000();


    ASSERT_EQ(rf == rf2, true);


    rf = ReferenceFrame::createBodyFixedSpice(1);
    ASSERT_EQ(rf.getType(), astro::ReferenceFrameType::BodyFixedNonRotating);
    rf = ReferenceFrame::createBodyFixedSpice(3);
    ASSERT_EQ(rf.getType(), astro::ReferenceFrameType::BodyFixedNonRotating);
    rf = ReferenceFrame::createBodyFixedSpice(6);
    ASSERT_EQ(rf.getType(), astro::ReferenceFrameType::BodyFixedNonRotating);
    rf = ReferenceFrame::createBodyFixedSpice(9);
    ASSERT_EQ(rf.getType(), astro::ReferenceFrameType::BodyFixedNonRotating);


    rf = ReferenceFrame::createBodyFixedSpice(301);
    ASSERT_EQ(rf.getType(), astro::ReferenceFrameType::BodyFixedRotating);
    rf = ReferenceFrame::createBodyFixedSpice(399);
    ASSERT_EQ(rf.getType(), astro::ReferenceFrameType::BodyFixedRotating);
    rf = ReferenceFrame::createBodyFixedSpice(699);
    ASSERT_EQ(rf.getType(), astro::ReferenceFrameType::BodyFixedRotating);





}

TEST_F(ReferenceFrameTest, OperatorTest)
{
    // CHeck that two J2000 frames gives a relative orientaion that equals identity
    astro::ReferenceFrame rf = ReferenceFrame::createJ2000();
    astro::ReferenceFrame rf2 = ReferenceFrame::createJ2000();
    astro::ReferenceFrame rf3 = ReferenceFrame::createBodyFixedSpice(399);
    astro::ReferenceFrame rf4 = ReferenceFrame::createBodyFixedSpice(399);



    ASSERT_EQ(rf, rf2);
    ASSERT_EQ(rf==rf3, false);
    ASSERT_EQ(rf3, rf4);


}


TEST_F(ReferenceFrameTest, BodyFixedFrameTest1)
{
    // OK for barycenters without loading any kernels
    // Also ok for many of the default frames wihtout loading kernels
    std::vector<int> ids = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 399, 301, 499, 401, 402,
        599, 501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516,
        699, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614,
            615, 616, 618, 
        799, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 
        899, 801, 802, 803, 804, 805, 806, 807, 808,
        999, 901};
    for(auto i : ids) {
        astro::ReferenceFrame rf2 = ReferenceFrame::createBodyFixedSpice(i);
        //std::cout << "ID: " << i << ", Frame " << rf2.getId() << ": " << rf2.getName() << std::endl;
        
        
        // This is some unrelated code to generate a list of planets for the settings files:
        std::string name = rf2.getName().substr(5);
        std::transform(name.begin(), name.end(), name.begin(), ::tolower);
        name.insert(0, 1, rf2.getName().substr(4,1)[0]);
        std::cout << "<celestialObject name=\"" << name << "\" id=\"" << i << "\" parentid=\"" << i/100 << "\"/>" << std::endl;
    }
}

TEST_F(ReferenceFrameTest, BodyFixedRotationTest1)
{
    // Removed this test, as the kernel was loade by some other test when runnign all tests
    // Calling getRotaion without any kernels loaded should throw exception:
    astro::ReferenceFrame rf_earth = ReferenceFrame::createBodyFixedSpice(399);
    mork::mat3d tip;
    //ASSERT_THROW(tip = rf_earth.getRotation(EphemerisTime(0)),
    //       astro::SpiceException);

    // This should be ok:
    astro::Spice().loadKernel("../data/spice/pck/pck00010.tpc");
    ASSERT_NO_THROW(tip = rf_earth.getRotationToJ2000(EphemerisTime(0))); 

    std::cout << tip << std::endl; 
} 

TEST_F(ReferenceFrameTest, FromStringTest1)
{
    astro::ReferenceFrame rf_J2000 = astro::ReferenceFrame::fromString("J2000");
    ASSERT_TRUE(rf_J2000.isJ2000());

    ASSERT_THROW(astro::ReferenceFrame err = astro::ReferenceFrame::fromString("TULL"), std::runtime_error);
}
