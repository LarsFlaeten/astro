#ifndef _ASTRO_RK1_4_H_
#define _ASTRO_RK1_4_H_

#include <vector>
#include "Time.h"
#include "State.h"
#include "ODE.h"
namespace astro
{

template<int N>
class RK
{
public:
    struct Result
    {
        State s;
        EphemerisTime et;
    };


    static Result doStep(const ODE& ode, const State& s, const EphemerisTime& et, const TimeDelta& dt);            

    static std::vector<Result> doSteps(const ODE& ode, const State& s, const EphemerisTime& et0, const EphemerisTime& et1, const TimeDelta& dt);

private:
    static const int n_stages;

};

template<int N>
const int RK<N>::n_stages=N;

template<int N>
typename RK<N>::Result RK<N>::doStep(const ODE& ode, const State& s, const EphemerisTime& et, const TimeDelta& dt)
{
    std::vector<double> a;
    std::vector< std::vector<double> > b;
    std::vector<double> c;
    switch (n_stages)
    {
    case 1:
        a = {0};
        b = { {0} };
        c = {1};
        break;
    case 2:
        a = {0, 1};
        b = { {0},
              {1}};
        c = {0.5, 0.5};
        break;

    case 3:
        a = {0, 0.5, 1};
        b = {   {0,   0},
                {0.5, 0.0},
                {-1.0, 2.0} };
        c = {1.0/6.0, 2.0/3.0, 1.0/6.0};
        break;

    case 4:
        a = {0, 0.5, 0.5, 1};
        b = { { 0.0,  0.0,  0.0 },
              { 0.5,  0.0,  0.0 },
              { 0.0,  0.5,  0.0 },
              { 0.0,  0.0,  1.0 }};
        c = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
        break;

    default:
        throw AstroException("RK Order needs to be <1..4>");
        break;    
    }



    double h = dt.value;
    double t_inner, ti = et.getETValue();
    astro::State s_inner, si = s;
    std::vector<astro::State> f;

    // Evaluate the time derivates at N points within the interval dt
    for(int i = 1; i <= n_stages; ++i)
    {
        t_inner = ti + a[i-1]*h;
        s_inner = si;
        for(int j = 1; j <= i-1; ++j)
        {   
            s_inner.r += (f[j-1].r)*(h*b[i-1][j-1]);
            s_inner.v += (f[j-1].v)*(h*b[i-1][j-1]);
        }
        f.push_back(ode.rates(EphemerisTime(t_inner), s_inner));
    }


    // prepare result:
    Result res;
    res.et = et + dt;
    res.s = si;
    for(int i = 0; i < 6; ++i)
    {
        res.s.r += f[i].r*(h*c[i]);
        res.s.v += f[i].v*(h*c[i]);
    }
    
	return res;

}
template<int N>
std::vector<typename RK<N>::Result> RK<N>::doSteps(const ODE& ode, const State& s, const EphemerisTime& et0, const EphemerisTime& et1, const TimeDelta& dt)
{
    std::vector<Result> res;
    // Initial value:
    res.push_back( RK<N>::Result({s, et0}));


    TimeDelta dti = dt;

    while(res.back().et < et1)
    {
        res.push_back(doStep(ode, res.back().s, res.back().et, dti));
        
        if(res.back().et + dti > et1)
            dti = et1 - res.back().et;
         
    }


    return std::move(res);
}            




}


#endif
