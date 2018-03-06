// Reference from some of the tests:
// [1]  Orbital Mechanics for Engineering Students, 2nd Edition, Howard D. Curtis


#include "../astro/Util.cpp"
#include <gtest/gtest.h>

#include <cmath>

class UtilTest : public ::testing::Test {

protected:
    UtilTest();

    virtual ~UtilTest();

    // Code here will be called immediately after the constructor (right
    // before each test).
    virtual void SetUp();

    // Code here will be called immediately after each test (right
    // before the destructor).
    virtual void TearDown();
};



UtilTest::UtilTest()
{

}

UtilTest::~UtilTest()
{

}

void UtilTest::SetUp()
{
}

void UtilTest::TearDown()
{
}

TEST_F(UtilTest, PiTest)
{
    double pi = astro::PI;
    ASSERT_EQ(pi, acos(-1.0));

    pi = astro::pi();
    ASSERT_EQ(pi, acos(-1.0));

    pi = astro::PI;
    ASSERT_EQ(pi, acos(-1.0));


}

TEST_F(UtilTest, DegreesRadiansTest)
{
    double deg, rad;
	
    deg = 0;
    rad = astro::radiansPerDegree()*deg;
    ASSERT_EQ(rad, 0.0);
    rad = astro::RADPERDEG*deg;
    ASSERT_EQ(rad, 0.0);

    deg = 90;
    rad = astro::radiansPerDegree()*deg;
    ASSERT_EQ(rad, astro::PI / 2.0);
    rad = astro::RADPERDEG*deg;
    ASSERT_EQ(rad, astro::PI / 2.0);


    for( auto i = 0; i <= 360; i++)
    {
        deg = static_cast<double>(i);
      
        rad = astro::RADPERDEG*deg;
        ASSERT_LT(fabs(rad - deg * astro::PI / 180.0), 1.0E-10);


        rad = deg * astro::PI / 180.0;
        deg = astro::DEGPERRAD*rad;
        ASSERT_LT(fabs(deg - static_cast<double>(i)), 1.0E-10);

    }
	
}

TEST_F(UtilTest, RADECTest)
{
    double dpr = astro::DEGPERRAD;
    // Ref [1], Example 4.1, plus two modifications to get negative dec, and rad in the other quadrant
    // (1)
    ork::vec3d r(-5368.0, -1784.0, 3691.0); // [km]
    double ra, dec, range;
    astro::vecToRaDec(r, &range, &ra, &dec); 
    
    ra = ra * dpr;
    dec = dec * dpr;
    ASSERT_LT(fabs(dec - 33.12), 1.0E-2);
    ASSERT_LT(fabs(ra - 198.4), 1.0E-1);

    // (2) - negative dec
    r = ork::vec3d(-5368.0, -1784.0, -3691.0);
    astro::vecToRaDec(r, &range, &ra, &dec);
    ra = ra * dpr;
    dec = dec * dpr;
	ASSERT_LT(fabs(-33.12 - dec), 1.0E-2);
    ASSERT_LT(fabs(ra - 198.4), 1.0E-1);

    // (3) - other quadrant of rad
    r = ork::vec3d(-5368.0, 1784.0, -3691.0);
    astro::vecToRaDec(r, &range, &ra, &dec);
    ra = ra * dpr;
    dec = dec * dpr;
    ASSERT_LT(fabs(-33.12 - dec), 1.0E-2);
    ASSERT_LT(fabs(ra - (360.0 - 198.4)), 1.0E-1);
}

TEST_F(UtilTest, NormalizeTest)
{
    double p2 = astro::TWOPI;

    double v = astro::wrap(0.0, 0.0, p2);
    ASSERT_EQ(v, 0.0);
    v = astro::wrap(p2-0.0000001, 0.0, p2);
    ASSERT_LT(fabs(v-p2), 0.000001);

    v = astro::wrap(3.1, 0.0, p2);
    ASSERT_EQ(v, 3.1);
    v = astro::wrap(2.7, 0.0, p2);
    ASSERT_EQ(v, 2.7);
    v = astro::wrap(1.44444, 0.0, p2);
    ASSERT_EQ(v, 1.44444);
    v = astro::wrap(6.11111, 0.0, p2);
    ASSERT_EQ(v, 6.11111);
    v = astro::wrap(5.99999, 0.0, p2);
    ASSERT_EQ(v, 5.99999);
    
    
    v = astro::wrap(45*p2+3, 0.0, p2);
    ASSERT_EQ(v, 3);
 
    v = astro::wrap(-45*p2+3, 0.0, p2);
    ASSERT_EQ(v, 3);
    
    v = astro::wrap(1000000000*p2+2.9999, 0.0, p2);
    ASSERT_LT(fabs(v- 2.9999), 1.0E-5);
 

}
