#include "RKF45.h"
#include "Exceptions.h"

#include <limits>
#include <algorithm>
#include <iostream>
#include <sstream>

namespace astro
{

const std::vector<double> RKF45::a = {0, 1./4., 3./8., 12./13., 1., 1./2.};

const std::vector< std::vector<double> > RKF45::b =
  { {     0,          0,          0,          0,         0, },
    {    1./4.,         0,          0,          0,         0 },
    {    3./32.,       9./32.,        0,          0,         0 },
    { 1932./2197., -7200./2197.,  7296./2197.,      0,         0 },
    {  439./216.,      -8.,      3680./513.,   -845./4104.,     0 },
    {   -8./27.,        2.,     -3544./2565.,  1859./4104.,  -11./40.}};
const std::vector<double> RKF45::c4 = {25./216., 0, 1408./2565., 2197./4104., -1./5., 0};
const std::vector<double> RKF45::c5 = {16./135., 0, 6656./12825., 28561./56430., -9./50., 2./55.}; 

const double RKF45::eps = std::numeric_limits<double>::epsilon();

const double RKF45::h_min = 16.0*RKF45::eps;

double RKF45::tol = 1.0E-8;

void    RKF45::setTolerance(double _tol)
{
    if(_tol < 0.0)
        throw AstroException("Zero or negative tolerance not allowed for RK methods");
    RKF45::tol = _tol;
}

RKF45::Result RKF45::doStep(const ODE& ode, const PosState& s, const EphemerisTime& et, const TimeDelta& dt)
{
    double h = dt.value;
    double t_inner, ti = et.getETValue();
    astro::PosState s_inner, si = s;
    std::vector<astro::PosState> f;

    // Evaluate the time derivates at six points within the interval dt
    for(int i = 1; i <= 6; ++i)
    {
        t_inner = ti + a[i-1]*h;
        s_inner = si;
        for(int j = 1; j <= i-1; ++j)
        {   
            s_inner.r += f[j-1].r*(h*b[i-1][j-1]);
            s_inner.v += f[j-1].v*(h*b[i-1][j-1]);
        }
        f.push_back(ode.rates(EphemerisTime(t_inner), s_inner));
    }

    // Compute maximum truncation error:
    astro::PosState te;
    te.r = vec3d::ZERO;
    te.v = vec3d::ZERO;
    for(int i = 0; i < 6; ++i)
    {
        
        te.r += f[i].r*(h*(c4[i]-c5[i]));
        te.v += f[i].v*(h*(c4[i]-c5[i]));
        
    }
    std::vector<double> te_abs = {fabs(te.r.x), fabs(te.r.y), fabs(te.r.z), fabs(te.v.x), fabs(te.v.y), fabs(te.v.z)};
    std::vector<double>::iterator te_maxi = std::max_element(te_abs.begin(), te_abs.end());
    double te_max = *te_maxi;


    // Compute allowable truncation error:
    std::vector<double> y_abs = {fabs(s.r.x), fabs(s.r.y), fabs(s.r.z), fabs(s.v.x), fabs(s.v.y), fabs(s.v.z)};
    std::vector<double>::iterator y_maxi = std::max_element(y_abs.begin(), y_abs.end());
    double te_allowed = std::max(*y_maxi,1.0) * tol;

    // Compute the fractional change in step size:
    double delta = pow(te_allowed/(te_max + eps), 1.0/5.0);
    
    // Update the time step:
    // We add a reduction factor here, since we hade the issue with
    // every other step being thrown away since h_next gave too high error
    // in the next step when h_next is trending downwards (after apoapsis).
    // This value must be checked with other orbits as well
    double h_next = std::min(0.93*delta*h, 4*h);
    if(h_next < h_min)
    {   
        std::stringstream ss;
        ss << "RKF45: Next Delta T fell below its minimim allowable value " << h_min << " at time: " << et.getETValue();
        throw AstroException(ss.str());
    } else if(te_max > te_allowed)
    {
        std::cout << "Oops... next time step will be: " << h_next << std::endl;
        //std::stringstream ss;
        //ss << "RKF45: Truncation error " << te_max << " is above allowed (" << te_allowed << ") at time: " << et.getETValue();
        //throw AstroException(ss.str());
        std::cout << "Skipping this step..." << std::endl;
        return Result({s, et, TimeDelta(h_next)});
    }

    // prepare result:
    Result res;
    res.dt_next = TimeDelta(h_next);
    res.et = et + dt;
    res.s = si;
    for(int i = 0; i < 6; ++i)
    {
        res.s.r += f[i].r*(h*c5[i]);
        res.s.v += f[i].v*(h*c5[i]);
    }
    
	return res;

}
 
std::vector<RKF45::Result> RKF45::doSteps(const ODE& ode, const PosState& s, const EphemerisTime& et0, const EphemerisTime& et1, const TimeDelta& dt)
{
    std::vector<Result> res;
    // Initial value:
    res.push_back( RKF45::Result({s, et0, dt}));


    while(res.back().et < et1)
    {
        res.push_back(doStep(ode, res.back().s, res.back().et, res.back().dt_next));
        
        if(res.back().et + res.back().dt_next > et1)
            res.back().dt_next = et1 - res.back().et;
         
    }


    return std::move(res);
}            


}
