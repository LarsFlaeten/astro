#include "Time.h"
#include "SpiceCore.h"

#include <cspice/SpiceUsr.h>


#include <mutex>
#include <sstream>
#include <iostream>

namespace astro
{


EphemerisTime::EphemerisTime()
    : et(0.0)
{

}

EphemerisTime::EphemerisTime(double _et)
    : et(_et)
{

}

EphemerisTime EphemerisTime::fromString(const std::string& datetime)
{
    double et;
    
    {
        std::lock_guard<std::mutex> lock(astro::Spice().mutex());
        str2et_c(datetime.c_str(), &et);
    }
	astro::Spice().checkError();
    return EphemerisTime(et);
}

EphemerisTime EphemerisTime::fromJED(double jed)
{
    // no mutex lock, only access constant spice functions    
    
    // REF: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/j2000_c.html
    double et = (jed - j2000_c()) * spd_c();
    
    return EphemerisTime(et);
}

EphemerisTime EphemerisTime::fromJDUTC(double jd)
{
    std::stringstream oss;    
    oss.setf(std::ios_base::fixed);
    oss << jd << " JD";    
    
    return fromString(oss.str());
}



double	EphemerisTime::getETValue() const
{
	return et;
}

std::string EphemerisTime::toISOUTCString(int prec) const
{
    // Sanitize input
    if(prec < 0) prec = 0;
    if(prec > 20 ) prec = 20; // Assume not sensible
    char    str[24+prec];
    {
        std::lock_guard<std::mutex>  lock(astro::Spice().mutex());
        et2utc_c(et, "ISOC", prec, 24+prec, str);
    }
    astro::Spice().checkError();

    return std::string(str);
}

std::string EphemerisTime::toJDUTCString(int prec) const
{
    // Sanitize input
    if(prec < 0) prec = 0;
    if(prec > 20 ) prec = 20; // Assume not sensible
    char    str[24+prec];
    {
        std::lock_guard<std::mutex>  lock(astro::Spice().mutex());
        et2utc_c(et, "J", prec, 24+prec, str);
    }
    astro::Spice().checkError();

    return std::string(str);
}

double  EphemerisTime::toJED() const
{
    // No mutex lock here, only call constant spice functions
    return et / spd_c() + j2000_c();
}

double  EphemerisTime::toJDUTC() const
{
    // Need mutex lock here, since we have to get DeltaET (=ET - UTC)
    // from the spice system
    double deltaET;
    {
        std::lock_guard<std::mutex>  lock(astro::Spice().mutex());
        deltet_c(et, "ET", &deltaET);
    }
    astro::Spice().checkError();

    double utc = et - deltaET;

    return utc / spd_c() + j2000_c();
}







}
