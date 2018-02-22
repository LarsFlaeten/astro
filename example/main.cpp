#include <iostream>

#include <astro/SpiceCore.h>
#include <astro/Time.h>

int main(int argc, char **argv)
{

    // This example shows normal time access and manipulation

    // Load the lepseconds kernel
    astro::Spice().loadKernel("../data/spice/lsk/naif0012.tls");

    // The basic time keeping entity used in astro is Ephemeris time,
    // Same as Spice. This value is represented as seconds past the
    // J2000 epoch in the time system known as Barycentric Dynamical Time
    // (TDB).
    //
    // Some examples of initialization
    astro::EphemerisTime et1; // et = 0 represents the J2000 epoch
    
    std::cout << "Date/Time in ISO: " << et1.toISOUTCString() << ", in Julian Day(UTC): " << et1.toJDUTCString() << std::endl;
    
    et1 = astro::EphemerisTime::fromString("2018 February 22, 20:04:00 UTC"); // The time this line was written
    std::cout << "Date/Time in ISO: " << et1.toISOUTCString() << ", in Julian Day(UTC): " << et1.toJDUTCString() << std::endl;
    
    et1 = astro::EphemerisTime::fromString("2017-12-24T17:00:00.12"); // ISO format
    std::cout << "Date/Time in ISO: " << et1.toISOUTCString() << ", in Julian Day(UTC): " << et1.toJDUTCString() << std::endl;
    
    et1 = astro::EphemerisTime::fromString("2451515.2981 JD"); // From Julian date
    std::cout << "Date/Time in ISO: " << et1.toISOUTCString() << ", in Julian Day(UTC): " << et1.toJDUTCString() << std::endl;
 
    // Init with J2000, shall give 2000 JAN 01 12:00:00
    et1 = astro::EphemerisTime::fromJDUTC(2451545.0);
    std::cout << "Date/Time in ISO: " << et1.toISOUTCString() << ", in Julian Day(UTC): " << et1.toJDUTCString() << std::endl;
    


    return 0;
}
