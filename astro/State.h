#ifndef _ASTRO_STATE_H_
#define _ASTRO_STATE_H_

#include <sstream>

// Decided to use ork::math. Other users need to install ork for this to work.
// TODO: Fix dependency in CMakeLists so that this is clear for users
#ifndef ORK_API
#define ORK_API
#endif
#include <ork/math/vec3.h>

#include <boost/operators.hpp>
#include <boost/numeric/odeint.hpp>

using namespace ork;

namespace astro
{


// The state of an object in a given reference frame
class State
    : boost::addable1< State ,
      boost::addable2< State , double ,
      boost::subtractable1< State,
      boost::subtractable2< State, double,
      boost::multipliable2< State, double > > > > >
{
public:
    vec3d   r; // Position [km]
    vec3d   v; // Orbital Velocity [km/s]

    State()
        : r(vec3d::ZERO), v(vec3d::ZERO)
    {  }

    State(const vec3d& _r, const vec3d& _v)
        : r(_r), v(_v)
    {  }

    State(const double val)
        : r(vec3d(val, val, val)), v(vec3d(val, val, val))
    {  }

    State& operator+=( const State& o)
    {
        r += o.r;
        v += o.v;
        return *this;
    }

    State& operator-=( const State& o)
    {
        r -= o.r;
        v -= o.v;
        return *this;
    }

    State& operator*=( const double a)
    {
        r *= a;
        v *= a;
        return *this;
    }

};

// only required for steppers with error control
State operator/( const State &p1 , const State &p2 );
State abs( const State &p );



std::ostream& operator << (std::ostream& os, const astro::State& s);



}

// also only for steppers with error control
namespace boost { namespace numeric { namespace odeint {
template<>
struct vector_space_norm_inf< astro::State >
{
    typedef double result_type;
    double operator()( const astro::State &p ) const
    {
        using std::max;
        using std::abs;
        return  max ( max( max( abs( p.r.x ) , abs( p.r.y ) ) , abs( p.r.z ) ),
                      max( max( abs( p.v.x ) , abs( p.v.y ) ) , abs( p.v.z ) ) );
    }
};
} } }




#endif
