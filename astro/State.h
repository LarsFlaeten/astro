#ifndef _ASTRO_STATE_H_
#define _ASTRO_STATE_H_

#include <sstream>

#include "ReferenceFrame.h"
#include "Time.h"


// Decided to use mork::math. Other users need to install mork for this to wmork.
// TODO: Fix dependency in CMakeLists so that this is clear for users
#include <mork/math/vec3.h>
#include <mork/math/quat.h>

#include <boost/operators.hpp>
#include <boost/numeric/odeint.hpp>

using namespace mork;

namespace astro
{


// The translative state of an object in a given reference frame
class PosState
    : boost::addable1< PosState ,
      boost::addable2< PosState , double ,
      boost::subtractable1< PosState,
      boost::subtractable2< PosState, double,
      boost::multipliable2< PosState, double > > > > >
{
public:
    vec3d   r; // Position [km]
    vec3d   v; // Orbital Velocity [km/s]

    PosState()
        : r(vec3d::ZERO), v(vec3d::ZERO)
    {  }

    PosState(const vec3d& _r, const vec3d& _v)
        : r(_r), v(_v)
    {  }

    PosState(const double val)
        : r(vec3d(val, val, val)), v(vec3d(val, val, val))
    {  }

    PosState& operator+=( const PosState& o)
    {
        r += o.r;
        v += o.v;
        return *this;
    }

    PosState& operator-=( const PosState& o)
    {
        r -= o.r;
        v -= o.v;
        return *this;
    }

    PosState& operator*=( const double a)
    {
        r *= a;
        v *= a;
        return *this;
    }

    // Transforms this state between reference frames at the given ET
    PosState    transform(const ReferenceFrame& fromFr, const ReferenceFrame toFr, const EphemerisTime& et);

};

// only required for steppers with error control
PosState operator/( const PosState &p1 , const PosState &p2 );
PosState abs( const PosState &p );



std::ostream& operator << (std::ostream& os, const astro::PosState& s);

// The rotation state of an object in a given reference coordinate system
// Usually it is given in the global/inertial frame
class RotState
/*    : boost::addable1< RotState ,
      boost::addable2< RotState , double ,
      boost::subtractable1< RotState,
      boost::subtractable2< RotState, double,
      boost::multipliable2< RotState, double > > > > >
*/
{
public:
    quatd   q; // Rotation quaternion in the given frame
    vec3d   w; // angular velocities in the given frame [radians/s]

    RotState()
        : q(quatd::ONE), w(vec3d::ZERO)
    {  }

    RotState(const quatd& _q, const vec3d& _w)
        : q(_q), w(_w)
    {  }

/*    RotState(const double val)
        : q(val, val, val, val), w(vec3d(val, val, val))
    {  }


    // Is this meaningful? We cant really add quats. But lets implement
    // Since the integrators need it..
    RotState& operator+=( const RotState& o)
    {
        q.x += o.q.x;
        q.y += o.q.y;
        q.z += o.q.z;
        q.w += o.q.w;
        w += o.w;
        return *this;
    }

    RotState& operator-=( const RotState& o)
    {
        q.x -= o.q.x;
        q.y -= o.q.y;
        q.z -= o.q.z;
        q.w -= o.q.w;
        w -= o.w;
        return *this;
    }

    // Again, not really meaningful, but used by the steppers
    RotState& operator*=( const double a)
    {
        q.x *= a;
        q.y *= a;
        q.z *= a;
        q.w *= a;
        w *= a;
        return *this;
    }
*/
};


std::ostream& operator << (std::ostream& os, const astro::RotState& s);



}

// Required by ODEINT steppers with error control
namespace boost { namespace numeric { namespace odeint {
template<>
struct vector_space_norm_inf< astro::PosState >
{
    typedef double result_type;
    double operator()( const astro::PosState &p ) const
    {
        using std::max;
        using std::abs;
        return  max ( max( max( abs( p.r.x ) , abs( p.r.y ) ) , abs( p.r.z ) ),
                      max( max( abs( p.v.x ) , abs( p.v.y ) ) , abs( p.v.z ) ) );
    }
};
} } }




#endif
