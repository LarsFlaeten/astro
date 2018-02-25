#include "Orbit.h"
#include "Util.h"
#include "SpiceCore.h"

#include <iostream>
#include <mutex>

#include <cspice/SpiceUsr.h>

#ifndef ORK_API
#define ORK_API
#endif
#include <ork/math/vec3.h>
#include <ork/math/mat3.h>
#include <ork/math/mat4.h>

using ork::vec3d;
using ork::mat3d;
using ork::mat4d;

namespace astro
{

OrbitElements OrbitElements::fromStateVectorOE(const State& state, double mu)
{
    if(mu < 0)
        throw AstroException("Non-positive mass give for primary body");

    // Method is from [1]
    OrbitElements oe;

    // 1. Calculate the distance:
    double r = state.r.length();
    if( r < 1.0E-5)
         throw AstroException("Radius of state vector near zero - Degenerate case");


    // 2. Scalar speed:
    double v = state.v.length();
    
    if( v < 1.0E-5)
        throw AstroException("Velocity near zero - Degenerate case");

    // 3. Radial velocity:
    // Note that:
    // - if v_r > 0, the sattelite is flying away from perigee
    // - if v_r < 0, it is flying towards perigee
    // v_r = R dot V/r
    double v_r = state.r.dotproduct(state.v/r);

    // 4. Specific angular momentum: H = R cross V
    vec3d H = state.r.crossProduct(state.v);

    // 5. Magnitude of angular momentum, the first orbital element:
    oe.h = H.length();
    if( oe.h < 1.0E-5)
        throw AstroException("Angular momentum near zero - Degenerate case");

    // 6. Inclination, the second orbital element
    // i = acos(h_z / h)
    // 0 < i < 90 is a prograde orbit
    // 90 < i < 180 is a retrograde orbit
    oe.i = acos(H.z / oe.h);

    // 7. Node line:
    // N = K cross H, K is the Z-unit vector
    vec3d N = vec3d(-H.y, H.x, 0.0);

    // 8. Magnitude of node line:
    double n = N.length();

    // 9. RA of the ascending node, the third orbital element
    //         { acos(N_x / n)          if N_y >= 0
    // Omega = { 
    //         { 360 - acos(N_x / n)    if N_y < 0
    oe.omega = acos(N.x / n);
    if(N.y < 0.0)
        oe.omega = 2.0*astro::PI - oe.omega;

    // 10. Eccentricity vector:
    vec3d R_prime = state.r*(v*v - mu/r);
    vec3d V_prime = state.v*(r*v_r);
    vec3d E = R_prime - V_prime;
    E *= (1.0 / mu); 

    // 11. Eccentricity, the fourth orbital element:
    // TODO: Check if eq. 4.11 [1] give improvements in speed. No, we need E for the argument of the perigee below..
    oe.e = E.length();

    // 12. Argument of perigee, fifth orbital element:
    //     { acos (N dot E / (n*e) )        if E_z >= 0
    // w = {
    //     { 360 - acos (N dot E / (n*e) )  if E_z <  0
    vec3d N_norm = N.normalize();
    vec3d E_norm = E.normalize();
    oe.w = acos(N_norm.dotproduct(E_norm));
    if(E.z < 0.0)
        oe.w = 2.0*astro::PI - oe.w;

    // 13. True anomaly, the sixth orbital element
    // theta = acos (E_norm dot R_norm)
    // theta = theta        if v_r >= 0
    // theta = 360 - theta  if v_r <  0
    vec3d R_norm = state.r.normalize();
    oe.theta = acos(E_norm.dotproduct(R_norm));
    if(v_r < 0.0)
        oe.theta = 2.0*astro::PI - oe.theta;

    // Additional elements:
    // Eccentric anomaly:
    double ea = 2.0*atan( tan(oe.theta / 2.0) / sqrt( (1 + oe.e) / (1 - oe.e) ) );

    // Mean anomaly:
    oe.M0 = ea - oe.e*sin(ea);

    // Periapsis distance:
    oe.rp = oe.h*oe.h / mu * (1.0 / (1.0 + oe.e));

    return oe;
}

#if 0
OrbitElements OrbitElements::fromstateVectorOEOpt(const State& state, double mu)
{
    // Method is from [1]
    OrbitElements oe;

    // 1. Calculate the distance:
    double r = state.r.length();

    // 2. Scalar speed:
    double v = state.v.length();

    // 3. Radial velocity:
    // Note that:
    // - if v_r > 0, the sattelite is flying away from perigee
    // - if v_r < 0, it is flying towards perigee
    // v_r = R dot V/r
    double v_r = state.r.dotproduct(state.v/r);

    // 4. Specific angular momentum: H = R cross V
    vec3d H = state.r.crossProduct(state.v);

    // 5. Magnitude of angular momentum, the first orbital element:
    oe.h = H.length();

    // 6. Inclination, the second orbital element
    // i = acos(h_z / h)
    // 0 < i < 90 is a prograde orbit
    // 90 < i < 180 is a retrograde orbit
    oe.i = acos(H.z / oe.h);

    // 7. Node line:
    // N = K cross H, K is the Z-unit vector
    //vec3d N = vec3d::UNIT_Z.crossProduct(H);
    //vec3d N = vec3d(-H.y, H.x, 0.0);

    // 8. Magnitude of node line:
    //double n = N.length();
    double n = sqrt(H.y*H.y + H.x*H.x);

    // 9. RA of the ascending node, the third orbital element
    //         { acos(N_x / n)          if N_y >= 0
    // Omega = { 
    //         { 360 - acos(N_x / n)    if N_y < 0
    //oe.omega = acos(N.x / n);
    oe.omega = acos(-H.y / n );

    //if(N.y < 0.0)
    if(H.x < 0.0)
        oe.omega = 2.0*astro::PI - oe.omega;

    // 10. Eccentricity vector:
    vec3d E = state.r*(v*v - mu/r) - state.v*(r*v_r);
    E *= (1.0 / mu); 

    // 11. Eccentricity, the fourth orbital element:
    // TODO: Check if eq. 4.11 [1] give improvements in speed. No, we need E for the argument of the perigee below..
    oe.e = E.length();

    // 12. Argument of perigee, fifth orbital element:
    //     { acos (N dot E / (n*e) )        if E_z >= 0
    // w = {
    //     { 360 - acos (N dot E / (n*e) )  if E_z <  0
    //vec3d N_norm = N.normalize();
    //N /= n;
    //vec3d E_norm = E.normalize();
    E /= oe.e; // Don't need E in its original form anymore
    //oe.w = acos(N_norm.dotproduct(E_norm));
    //oe.w = acos(N.dotproduct(E));
    oe.w = acos(-H.y*E.x/n + H.x*E.y/n);
    if(E.z < 0.0)
        oe.w = 2.0*astro::PI - oe.w;

    // 13. True anomaly, the sixth orbital element
    // theta = acos (E_norm dot R_norm)
    // theta = theta        if v_r >= 0
    // theta = 360 - theta  if v_r <  0
    vec3d R_norm = state.r.normalize();
    oe.theta = acos(E.dotproduct(R_norm));
    if(v_r < 0.0)
        oe.theta = 2.0*astro::PI - oe.theta;


    return oe;
}
#endif

State   OrbitElements::toStateVectorOE(double mu)
{
    // From [1], Algorithm 4.5:

    // Calculate Rxp and Vxp in perifocal frame:
    double costheta = cos(theta);
    double sintheta = sin(theta);
    vec3d Rxp(costheta, sintheta, 0.0);
    Rxp *= (h * h)/mu * (1 / ( 1 + e*costheta ));
    //std::cout << "RXP: [" << Rxp.x << ", " << Rxp.y << ", " << Rxp.z << "]" << std::endl;


    vec3d Vxp(-sintheta, e + costheta, 0.0);
    Vxp *= mu / h;
    //std::cout << "VXP: [" << Vxp.x << ", " << Vxp.y << ", " << Vxp.z << "]" << std::endl;


    // Construct Rotatation matrix from perifocal to geocentric equatorial:
    // MOdified with according to [2] as [¡] is inconsistent with rotations
    mat3d R3_w = mat4d::rotatez(-w*DEGPERRAD).mat3x3();
    mat3d R1_i = mat4d::rotatex(-i*DEGPERRAD).mat3x3();
    mat3d R3_Omega = mat4d::rotatez(-omega*DEGPERRAD).mat3x3();
    
    // Had to use the transpose rotation matrices from ork, since [1] uses
    // the transpose equivalent of these. Why do this?? The ork ones seem right, they match
    // With rotatino matrices online for RH coordinate systems.. It doesnt seem right what [1] does..
    // But at least all the tests pass by using the transposes.
    // TODO: Look into this later...
    //mat3d Q_X_to_xp = R3_w * R1_i * R3_Omega;
    mat3d Q_X_to_xp = R3_w * R1_i * R3_Omega;
    mat3d Q_xp_to_X = Q_X_to_xp.transpose();

    State state;

    state.r = Q_xp_to_X * Rxp;
    state.v = Q_xp_to_X * Vxp;


    return state;
}

OrbitElements OrbitElements::fromStateVectorSpice(const State& state, double mu)
{
    astro::Spice(); // Init spice
    
    double et = 0.0;
    double elts[8];

    // need mutex?
    // This function shold not look up any internal states
    // but uses Spice error handling which might not be thread safe.
    // Did some checks, this method is 5-6 times slower than the OE counterpart. is it checkError
    // Who is response for this, not the OSCELT algorithm.
    // TODO: Optimize checkError...
    
    {
        std::lock_guard<std::mutex> lock(astro::Spice().mutex());
        oscelt_c(&(state.r.x), et, mu, elts);
    }
    astro::Spice().checkError();


    OrbitElements oe;
    oe.rp = elts[0];
    oe.e = elts[1];
    oe.i = elts[2];
    oe.omega = elts[3];
    oe.w = elts[4];

    // angular momentum:
    oe.h = state.r.crossProduct(state.v).length();

    // Mean anomaly 
    // [1], (4.13b)
    double r = state.r.length();
    double v_r  = state.r.dotproduct(state.v) / r;
    oe.theta = acos( (1.0 / oe.e) * (oe.h*oe.h/(mu*r ) - 1) );
    if(v_r < 0)
        oe.theta = 2.0*astro::PI - oe.theta;
    oe.M0 = elts[5];

    return oe;
}


State   OrbitElements::toStateVectorSpice(double mu)
{

}
}
