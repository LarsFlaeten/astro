#ifndef _ASTRO_TIME_H_
#define _ASTRO_TIME_H_

#include <string>

namespace astro
{

struct TimeDelta
{
    TimeDelta(double dt)
       : value(dt)
    {

    }

    double value;
};

class EphemerisTime
{
public:
    // Initializes with ET = 0 (J2000)
    EphemerisTime();

    // Initialize with given ET.
    EphemerisTime(double et); 

    // Normally this statuc function is used for initialization
    static EphemerisTime fromString(const std::string& timedate);
    
    // ..or this one, from julian day based on TDB/ET
    static EphemerisTime fromJED(double jed);

    // ..or this one, from UTC based JD
    static EphemerisTime fromJDUTC(double jd);


	// returns the numerical ET value
	double getETValue() const;

    // Returns a string representation of the date&time in UTC ISO format
    // prec i precision in the seconds part
    std::string toISOUTCString(int prec = 0) const;

    // Returns a string representation of the date&time in JD UTC format
    // prec is precision in the day part
    std::string toJDUTCString(int prec = 5) const;

    // Returns the julian day (based on TDB/ET or UTC) equivalent of this et:
    double      toJED() const;
    double      toJDUTC() const;

    bool operator==(const EphemerisTime& other)
    {
        return et == other.et ? true : false;
    }

    void operator+=(const TimeDelta dt)
    {
        et += dt.value;
    }
    void operator-=(const TimeDelta dt)
    {
        et -= dt.value;
    }


private:
    // The ephemeris time
    double et;

};




}


#endif
