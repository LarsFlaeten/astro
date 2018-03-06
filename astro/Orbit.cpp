#include "Orbit.h"
#include "Util.h"
#include "SpiceCore.h"

#include <iostream>
#include <sstream>
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

OrbitElements OrbitElements::fromStateVectorOE(const State& state, const EphemerisTime& epoch, double mu)
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
    if(n>0.0)
        oe.omega = acos(N.x / n);
    else
        oe.omega = 0.0;
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
    vec3d E_norm = E.normalize();
    if( n > 0.0)
    {
        vec3d N_norm = N.normalize();
       oe.w = acos(N_norm.dotproduct(E_norm));
    } else
        oe.w = 0.0;

    if(E.z < 0.0)
        oe.w = 2.0*astro::PI - oe.w;

    // 13. Mean anomaly from true anomaly, the sixth orbital element
    // theta = acos (E_norm dot R_norm)
    // theta = theta        if v_r >= 0
    // theta = 360 - theta  if v_r <  0
    vec3d R_norm = state.r.normalize();
    double theta = acos(E_norm.dotproduct(R_norm));
    if(v_r < 0.0)
        theta = 2.0*astro::PI - theta;

    // Mean anomaly
    oe.M0 = meanAnomalyFromTrueAnomaly(theta, oe.e);

    /// Assign epoch and mu
    oe.mu = mu;
    oe.epoch = epoch;
    
    // Compute derived, but constant paramaters:
    oe.computeDerivedQuantities();

    return oe;

}

void OrbitElements::computeDerivedQuantities()
{
    // Periapsis distance:
    rp = h*h / mu * (1.0 / (1.0 + e));
    
    // Apoapsis distance:
    ap = h*h / mu * (1.0 / (1.0 - e));

    // Semimajor axis
    a = 0.5*(rp + ap);

    // Period (a negative number if e >= 1)
    if( e < 1.0)
        T = astro::TWOPI / sqrt(mu)*pow(a, 3.0/2.0);
    else
        T = -1.0;

    // Mean motion
    n = sqrt(mu/pow(fabs(a), 3.0));
}

double OrbitElements::meanAnomalyFromTrueAnomaly(double trueAnomaly, double e)
{
    if( e < 0.0 )
        throw astro::AstroException("ERROR, negative eccentricity is not allowed");
	
    if (e < 1.0)
    {
        // Eccentric anomaly:
        double E = 2.0*atan( tan(trueAnomaly / 2.0) / sqrt( (1 + e) / (1 - e) ) );

        // Mean anomaly:
        return E - e*sin(E);
    } else
    {
        // http://control.asu.edu/Classes/MAE462/462Lecture05.pdf
        // Hyperbolic anomaly:
        double H = 2.0*atanh( (sqrt( (e-1) / (e + 1)  ) ) * tan(trueAnomaly/2.0) );
    
        return e*sinh(H) - H;

    }
}

double OrbitElements::trueAnomalyFromMeanAnomaly(double M, double e)
{
    if( e < 0.0 )
        throw astro::AstroException("ERROR, negative eccentricity is not allowed");
	
	std::pair<double,int> res;
    res = Kepler2(M, e);

	// calculate true anomaly from eccentric anomaly
	// https://en.wikipedia.org/wiki/True_anomaly
    if( e < 0.9998)
    {
        double E = res.first;
	    double cosE2 = cos(0.5*E);
        double sinE2 = sin(0.5*E);
	    double theta = 2.0*atan2(sqrt(1.0+e)*sinE2, sqrt(1.0-e)*cosE2);

	    return theta;
    }
    else
    {
        double H = res.first;
        double theta = 2.0*atan(sqrt(( e + 1.0 )/(e-1.0))*tanh(0.5*H));
        return theta;
    }
}
    


std::pair<double,int> OrbitElements::Kepler1(double M, double e)
{
    if ( e < 0.0 )
        throw astro::AstroException("ERROR, negative eccentricity is not allowed");
	if ( e >= 1.0 )
		throw astro::AstroException("ERROR, Kepler2 is for e < 1.0");


    // Reduce Anomaly to [0..2PI]
    double Mw = astro::wrap(M, 0.0, astro::TWOPI);

	double E0 = Mw; // First guess
    double En;
    int it = 0;

    while(true)
    {
        ++it;
        En = Mw + e*sin(E0);
        if(fabs(En-E0)<1.0E-10)
            break;
        if(it>= KEPLER_MAX_ITERATIONS)
            throw astro::AstroException("Kepler1 did not converge within max interations");
        E0 = En;
    }

    

    return std::make_pair(En,it);
}

std::pair<double,int> OrbitElements::Kepler2(double M, double e)
{
    if ( e < 0.0 )
        throw astro::AstroException("ERROR, negative eccentricity is not allowed");
	

    if ( e < 0.9998 )
	{
        // Reduce Anomaly to [0..2PI]
        double Mw = astro::wrap(M, 0.0, astro::TWOPI);

   	    double E0 = Mw; // First guess
        double En;
        int it = 0;

       
        while(true)
        {
            ++it;
        
            // Many references for this, e.g:
            // https://www.csun.edu/~hcmth017/master/node16.html
            En = E0 - (E0 - e*sin(E0) - Mw)/(1.0 - e*cos(E0));
        

            if(fabs(En-E0)<KEPLER_TOLERANCE)
                break;
            if(it>= KEPLER_MAX_ITERATIONS)
            {
                std::stringstream oss;
                oss << "Kepler2 did not converge within max interations(" << astro::KEPLER_MAX_ITERATIONS << "), e=" << e << ", M=" << Mw;
                
                throw astro::AstroException(oss.str());
            }
            
            E0 = En;
        }

        return std::make_pair(En,it);
    }
    else if (e > 1.0)
    {
        double H0 = M;
        double Hn;
        int it = 0; 
        while(true)
        {
            ++it;
            
            // http://control.asu.edu/Classes/MAE462/462Lecture05.pdf
            Hn = H0 + (M - e*sinh(H0) + H0)/(e*cosh(H0) - 1.0);
            

            if(fabs(Hn-H0)<KEPLER_TOLERANCE)
                break;
            if(it>= KEPLER_MAX_ITERATIONS)
            {
                std::stringstream oss;
                oss << "Kepler2 did not converge within max interations(" << astro::KEPLER_MAX_ITERATIONS << "), e=" << e << ", M=" << M;
                
                throw astro::AstroException(oss.str());
            }
            H0 = Hn;
        }   

        return std::make_pair(Hn,it);
 
    }
    else
     	throw astro::AstroException("ERROR, Kepler2 is for e < 0.9998 or e > 1.0");



}

State   OrbitElements::toStateVectorOE(const EphemerisTime& et)
{

    // Get the Mean Anomaly:
    double t0 = this->epoch.getETValue();
    double t = et.getETValue();
    double M = M0 + this->n*(t-t0);

    // True anomaly:
    double theta = trueAnomalyFromMeanAnomaly(M, this->e);

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

OrbitElements OrbitElements::fromStateVectorSpice(const State& state, const EphemerisTime& epoch, double mu)
{
    astro::Spice(); // Init spice
    
    double et = epoch.getETValue();
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

    // Mean anomaly at epoch 
    oe.M0 = elts[5];

    // Assign epoch and mu
    oe.mu = mu;
    oe.epoch = epoch;

    // Compute derived quantities:
    oe.computeDerivedQuantities();


    return oe;
}


State   OrbitElements::toStateVectorSpice(const EphemerisTime& et)
{
    double elts[8];

    elts[0] = rp;
    elts[1] = e;
    elts[2] = i;
    elts[3] = omega;
    elts[4] = w;
    elts[5] = M0;
    elts[6] = epoch.getETValue();
    elts[7] = mu;

    //for(auto i = 0; i < 8; ++i)
    //    std::cout << elts[i] << std::endl;

    double st[6];
    {
        std::lock_guard<std::mutex> lock(astro::Spice().mutex());
        conics_c(elts, et.getETValue(), st);
    }
    astro::Spice().checkError();

    astro::State state;
    state.r.x = st[0];
    state.r.y = st[1];
    state.r.z = st[2];
    state.v.x = st[3];
    state.v.y = st[4];
    state.v.z = st[5];
    
    return state;
}
}
