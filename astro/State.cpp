#include "State.h"

namespace astro
{

std::ostream& operator << (std::ostream& os, const astro::PosState& s)
{
    os << "r: (" << s.r.x <<", " << s.r.y << ", " << s.r.z << ") [km]\n";
    os << "v: (" << s.v.x <<", " << s.v.y << ", " << s.v.z << ") [km/s]\n";
    return os;
}   


PosState operator/(const PosState& p1, const PosState& p2)
{
    return PosState( p1.r/p2.r , p1.r/p2.r);
}


PosState abs(const PosState& p)
{
    return PosState( vec3d(std::abs(p.r.x) , std::abs(p.r.y) , std::abs(p.r.z) ),
					vec3d(std::abs(p.v.x) , std::abs(p.v.y) , std::abs(p.v.z)));

}



}
