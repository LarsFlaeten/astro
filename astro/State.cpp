#include "State.h"

namespace astro
{

std::ostream& operator << (std::ostream& os, const astro::State& s)
{
    os << "r: (" << s.r.x <<", " << s.r.y << ", " << s.r.z << ") [km]\n";
    os << "v: (" << s.v.x <<", " << s.v.y << ", " << s.v.z << ") [km/s]\n";
    return os;
}   

std::ostream& operator << (std::ostream& os, const astro::StateDot& s)
{
    os << "v: (" << s.v.x <<", " << s.v.y << ", " << s.v.z << ") [km/s]\n";
    os << "a: (" << s.a.x <<", " << s.a.y << ", " << s.a.z << ") [km²/s]\n";
    return os;
}   




}
