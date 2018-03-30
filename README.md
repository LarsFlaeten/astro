# astro - An astrodynamics toolbox

WORK IN PROGRESS!

Astro is a collection of astrodynamics-related classes and helper functions. In addition it wraps some functionality of NAIF's SPICE toolkit, but provides an interface in modern C++ fashion.

Math library is based on the Ork math templates.

Put cspice headers and libcspice.a in `libraries/cspice/` (not included in this repo)

An example app shows some of the functionality. The test cases also shows a lot of the functionality of the toolbox.

# Examples

## Example 1 - Time

The basic time keeping entity used in astro is Ephemeris time, which is the same a in the NAIF Spice toolkit. This value is represented as seconds past the J2000 epoch in the time system known as Barycentric Dynamical Time (BDT). In order to convert to UTC and other more human readable time formats, we use Spice to keep track of leap seconds. The leap seconds data is loaded into Spice / astro in the beginning of an application if time conversion is needed:

```
// Load the lepseconds kernel
astro::Spice().loadKernel("../data/spice/lsk/naif0012.tls");
```

Some examples of Ephemeris Time initialization and conversion:

```
astro::EphemerisTime et1; // et = 0 is default and represents the J2000 epoch
std::cout << "Date/Time in ISO: " << et1.toISOUTCString() << ", in Julian Day(UTC): " << et1.toJDUTCString() << std::endl;
```
Output:
```
Date/Time in ISO: 2000-01-01T11:58:56, in Julian Day(UTC): JD 2451544.99926
```

Initialization from UTC string:

```
et1 = astro::EphemerisTime::fromString("2018 February 22, 20:04:00 UTC");
std::cout << "Date/Time in ISO: " << et1.toISOUTCString() << ", in Julian Day(UTC): " << et1.toJDUTCString() << std::endl;
```
Output:
```
Date/Time in ISO: 2018-02-22T20:04:00, in Julian Day(UTC): JD 2458172.33611
```

Initialization from ISO string:


```
et1 = astro::EphemerisTime::fromString("2017-12-24T17:00:00.12"); // ISO format
std::cout << "Date/Time in ISO: " << et1.toISOUTCString() << ", in Julian Day(UTC): " << et1.toJDUTCString() << std::endl;
```
Output:
```
Date/Time in ISO: 2017-12-24T17:00:00, in Julian Day(UTC): JD 2458112.20833
```

Initialization from string with Julian Day:


```
et1 = astro::EphemerisTime::fromString("2451515.2981 JD"); // From Julian date
std::cout << "Date/Time in ISO: " << et1.toISOUTCString() << ", in Julian Day(UTC): " << et1.toJDUTCString() << std::endl;
```
Output:
```
Date/Time in ISO: 1999-12-02T19:09:16, in Julian Day(UTC): JD 2451515.29810
```

Initialization with Julian Day:
```
// Init with J2000, shall give 2000 JAN 01 12:00:00
et1 = astro::EphemerisTime::fromJDUTC(2451545.0);
std::cout << "Date/Time in ISO: " << et1.toISOUTCString() << ", in Julian Day(UTC): " << et1.toJDUTCString() << std::endl;
```
Output:
```
Date/Time in ISO: 2000-01-01T12:00:00, in Julian Day(UTC): JD 2451545.00000
```

## Example 2 - Orbital elements and simple orbits

This example shows how to establish a Keplerian (unperturbed, 2-body) orbit from a state vector.

```
astro::PosState   state;
state.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
state.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]

astro::EphemerisTime et(0); // et = 0 represents the J2000 epoch
double mu_earth = 398600.0;

// Convert to Keplerian orbital elements for this epoch
astro::OrbitElements oe = astro::OrbitElements::fromStateVector(state, et, mu_earth);

// Generate a simple orbit from these elements:
astro::SimpleOrbit orbit1(oe);

// Dump orbit elemements to stdout:
std::cout << orbit1.getOrbitElements() << std::endl;
```
Output:
```
Angular momentum:    58311.7 [km²/s]
Inclination:         153.249 [Deg]
RA of the asc. node: 255.279 [Deg]
Eccentricity:        0.171212 [-]
Argument of perigee: 20.0683 [Deg]
Mean anomaly @epoch: 20.0709 [Deg]
Epoch:               0 [seconds]
mu:                  398600 [km³/s²]
Periapsis distance:  7283.46 [km]
Apoapsis distance:   10292.7 [km]
Semimajor axis:      8788.1 [km]
Period:              8198.86 [s]
Mean motion:         0.000766349 [rad/s]
```

## Example 3 - Orbit positions as a function of time:

This example shows how to establish a keplerian orbit from a state vector, and how to retrieve orbit positions as a function of time.
```
astro::PosState   state;
state.r = vec3d(-6045.0, -3490.0, 2500.0);  //[km]
state.v = vec3d(-3.457, 6.618, 2.533);      //[km/s]

astro::EphemerisTime et(0); // et = 0 represents the J2000 epoch
double mu_earth = 398600.0;

// Convert to Keplerian orbital elements for this epoch
astro::OrbitElements oe = astro::OrbitElements::fromStateVector(state, et, mu_earth);

// Generate a simple orbit from these elements:
astro::SimpleOrbit orbit1(oe);

astro::TimeDelta dt(orbit1.getPeriod()/40.0);

for(int i = 0; i <= 40; ++i)
{
    // Get the state and print selected elements (position)
    auto state = orbit1.getState(et);
    std::cout << et.getETValue() << "\t" << state.r.x << "\t" << state.r.y << "\t" << state.r.z << std::endl;

    // Advance EphemerisTime with 1/40th period:
    et += dt;

}
```
The orbit is shown in the follofing plot:

![alt text](https://raw.githubusercontent.com/LarsFlaeten/astro/master/web/example3.png "Orbit from example 3 plotted aronud primary (Earth)")

