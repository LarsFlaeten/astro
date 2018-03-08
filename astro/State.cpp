#include "State.h"

std::ostream& operator << (std::ostream& os, const astro::State& s)
{
    os << "R: (" << s.r.x <<", " << s.r.y << ", " << s.r.z << ") [km]\n";
    os << "V: (" << s.v.x <<", " << s.v.y << ", " << s.v.z << ") [km/s]\n";
    return os;
}   

void print(const astro::State& s)
{
    std::cout << "R: (" << s.r.x <<", " << s.r.y << ", " << s.r.z << ") [km]\n";
    std::cout << "V: (" << s.v.x <<", " << s.v.y << ", " << s.v.z << ") [km/s]\n";
    std::cout << std::endl;
}


