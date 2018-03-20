#include "../astro/Time.cpp"
#include <gtest/gtest.h>

#include <cmath>

class TimeTest : public ::testing::Test {

protected:
    TimeTest();

    virtual ~TimeTest();

    // Code here will be called immediately after the constructor (right
    // before each test).
    virtual void SetUp();

    // Code here will be called immediately after each test (right
    // before the destructor).
    virtual void TearDown();
};



TimeTest::TimeTest()
{

}

TimeTest::~TimeTest()
{

}

void TimeTest::SetUp()
{
    // Need the leapseconds kernel:
    astro::Spice().loadKernel("../data/spice/lsk/naif0012.tls");
}

void TimeTest::TearDown()
{
    // TODO: Remove the LS kernel again..
}

TEST_F(TimeTest, ETInitializationTest)
{

    astro::EphemerisTime et;
    ASSERT_EQ(et.getETValue(), 0.0);
    
    // Spice docs:
    //
    // ISO T:
    ASSERT_NO_THROW(et = astro::EphemerisTime::fromString("1996-12-18T12:28:28"));
    ASSERT_NO_THROW(et = astro::EphemerisTime::fromString("1986-01-18T12"));
    ASSERT_NO_THROW(et = astro::EphemerisTime::fromString("1986-01-18T12:19"));
    ASSERT_NO_THROW(et = astro::EphemerisTime::fromString("1986-01-18T12:19:52.18"));
    ASSERT_NO_THROW(et = astro::EphemerisTime::fromString("1995-08T18:28:12"));
    ASSERT_NO_THROW(et = astro::EphemerisTime::fromString("1995-18T"));
    // various strings that should be possible to parse
    ASSERT_NO_THROW(et = astro::EphemerisTime::fromString("Tue Aug  6 11:10:57  1996"));
   	ASSERT_NO_THROW(et = astro::EphemerisTime::fromString("1 DEC 1997 12:28:29.192"));
   	ASSERT_NO_THROW(et = astro::EphemerisTime::fromString("2/3/1996 17:18:12.002"));
   	ASSERT_NO_THROW(et = astro::EphemerisTime::fromString("Mar 2 12:18:17.287 1993"));
   	ASSERT_NO_THROW(et = astro::EphemerisTime::fromString("1992 11:18:28  3 Jul"));
   	ASSERT_NO_THROW(et = astro::EphemerisTime::fromString("June 12, 1989 01:21"));
   	ASSERT_NO_THROW(et = astro::EphemerisTime::fromString("1978/3/12 23:28:59.29"));
   	ASSERT_NO_THROW(et = astro::EphemerisTime::fromString("17JUN1982 18:28:28"));
   	ASSERT_NO_THROW(et = astro::EphemerisTime::fromString("13:28:28.128 1992 27 Jun"));
   	ASSERT_NO_THROW(et = astro::EphemerisTime::fromString("1972 27 jun 12:29"));
   	ASSERT_NO_THROW(et = astro::EphemerisTime::fromString("'93 Jan 23 12:29:47.289"));
   	ASSERT_NO_THROW(et = astro::EphemerisTime::fromString("27 Jan 3, 19:12:28.182"));
   	ASSERT_NO_THROW(et = astro::EphemerisTime::fromString("23 A.D. APR 4, 18:28:29.29"));
   	ASSERT_NO_THROW(et = astro::EphemerisTime::fromString("18 B.C. Jun 3, 12:29:28.291"));
   	ASSERT_NO_THROW(et = astro::EphemerisTime::fromString("29 Jun  30 12:29:29.298"));
   	ASSERT_NO_THROW(et = astro::EphemerisTime::fromString("29 Jun '30 12:29:29.298"));
    
    // Should not accept empty string
    ASSERT_THROW(et = astro::EphemerisTime::fromString(""), astro::SpiceException);

	
	
}

TEST_F(TimeTest, EqualityTest)
{
    ASSERT_EQ(astro::EphemerisTime(0.0)==astro::EphemerisTime(0.0), true);
    ASSERT_EQ(astro::EphemerisTime(0.0)==astro::EphemerisTime(1.0), false);
    ASSERT_EQ(astro::EphemerisTime(1.0)==astro::EphemerisTime(1.0), true);
    ASSERT_EQ(
        astro::EphemerisTime::fromString("29 Jun '30 12:29:29.298")==
        astro::EphemerisTime::fromString("29 Jun '30 12:29:29.298"),
        true);
    ASSERT_EQ(
        astro::EphemerisTime::fromString("29 Jun '30 12:29:29.298")==
        astro::EphemerisTime(0.0),
        false);





}

TEST_F(TimeTest, Operatorstest)
{
    astro::EphemerisTime et0(0.0);
    astro::EphemerisTime et1(1.0);
    astro::TimeDelta dt(1.0);

    ASSERT_EQ(et1 > et0, true);
    ASSERT_EQ(et1 < et0, false);
    ASSERT_EQ(et0 > et1, false);
    ASSERT_EQ(et0 < et1, true);


    ASSERT_EQ(et1==et0, false);
    
    astro::EphemerisTime etn = et0 + dt;
    ASSERT_EQ(etn, et1);

    et0 += dt;
    ASSERT_EQ(et0, et1);


    ASSERT_EQ(et0==et1, true);
    
    et0 -= dt;
    ASSERT_EQ((et1 - et0).value, dt.value);

}


TEST_F(TimeTest, ET2UTCTest)
{

        astro::EphemerisTime et1, et2;

        // JD:
        // Do conversion from string to et and back to string
        std::string jd = "JD 2446533.1883428";
        et1 = astro::EphemerisTime::fromString(jd);
        std::string jd2 = et1.toJDUTCString(7);
        ASSERT_STREQ(jd.c_str(), jd2.c_str());        
        // ISO
        // Do conversion from string to et and back to string
        std::string iso = "1987-04-12T16:31:12.814";
        et2 = astro::EphemerisTime::fromString(iso);
        std::string iso2 = et2.toISOUTCString(3);
        ASSERT_STREQ(iso.c_str(), iso2.c_str());

        // Check negative precision is corrected to prec=0:
        // Note: JD is allways reported back with a dot, even with prec=0
        // ISO drops the dot if prec=0.
        
        jd = "JD 2446533.";
        et1 = astro::EphemerisTime::fromString(jd);
        jd2 = et1.toJDUTCString(-42);
        ASSERT_STREQ(jd.c_str(), jd2.c_str());

        iso = "1987-04-12T16:31:12";
        et2 = astro::EphemerisTime::fromString(iso);
        iso2 = et2.toISOUTCString(-655);
        ASSERT_STREQ(iso.c_str(), iso2.c_str());

}

TEST_F(TimeTest, JDConversions)
{
    // Init at ET=0 (J2000)
    astro::EphemerisTime et;
    double jed = et.toJED();
    ASSERT_EQ(jed, 2451545.0);
    

    // Check the UTC variants, these will be deltaET offset:
    double jd = et.toJDUTC();
    std::string jds = et.toJDUTCString(6);
    // jds and jd shall be allmost equal
    // Strip off the "JD ", and convert to double
    jds = jds.substr(3, jds.size());
    double jd2 = std::stod(jds);
    // difference should be less than 1.0E-6 (approx 0.08 seconds)
    ASSERT_LT(fabs(jd-jd2), 1.0E-6); 



}

TEST_F(TimeTest, TimeDeltaTest)
{
    astro::EphemerisTime et0(0);
    et0 += astro::TimeDelta(377.67876);
    ASSERT_EQ(et0.getETValue(),  377.67876);
    et0 -= astro::TimeDelta(377.67876);
    ASSERT_EQ(et0.getETValue(),  0.0);

    astro::EphemerisTime et1(9836474635);
    et1 += astro::TimeDelta(746.63546);
    ASSERT_EQ(et1.getETValue(), 9836475381.63546);
    et1 -= astro::TimeDelta(746.63546);
    ASSERT_EQ(et1.getETValue(),  9836474635);


 

 


}


