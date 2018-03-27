#include "State.h"

namespace astro
{

std::ostream& operator << (std::ostream& os, const astro::State& s)
{
    os << "r: (" << s.r.x <<", " << s.r.y << ", " << s.r.z << ") [km]\n";
    os << "v: (" << s.v.x <<", " << s.v.y << ", " << s.v.z << ") [km/s]\n";
    return os;
}   


State operator/(const State& p1, const State& p2)
{
    return State( p1.r/p2.r , p1.r/p2.r);
}


State abs(const State& p)
{
    return State( vec3d(std::abs(p.r.x) , std::abs(p.r.y) , std::abs(p.r.z) ),
					vec3d(std::abs(p.v.x) , std::abs(p.v.y) , std::abs(p.v.z)));

}



}
