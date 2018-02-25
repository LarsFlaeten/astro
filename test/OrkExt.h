#ifndef ORK_EXT
#define ORK_EXT
#include <sstream>


#ifndef ORK_API
#define ORK_API
#endif
#include <ork/math/vec3.h>
#include <ork/math/vec4.h>
#include <ork/math/mat3.h>
#include <ork/math/mat4.h>


namespace ork
{

// Printing a vec3<T> etc..
template < class T >
std::ostream& operator << (std::ostream& os, const ork::vec3<T>& v) 
{
        os << "[" << v.x << " " << v.y << " " << v.z << "]";
        return os;
}

template < class T >
std::ostream& operator << (std::ostream& os, const ork::vec4<T>& v) 
{
        os << "[" << v.x << " " << v.y << " " << v.z << " " << v.w << "]";
        return os;
}


template < class T >
std::ostream& operator << (std::ostream& os, const ork::mat4<T>& m) 
{
    ork::vec4<T> v;
    for(int i = 0; i < 4; i++)
    {
        v = m[i];
        os << "[" << v.x << " " << v.y << " " << v.z << " " << v.w << "]\n";
    }
    
    return os;
}

template < class T >
std::ostream& operator << (std::ostream& os, const ork::mat3<T>& m) 
{
    ork::vec3<T> v;
    for(int i = 0; i < 3; i++)
    {
        v = m[i];
        os << "[" << v.x << " " << v.y << " " << v.z << "]\n";
    }
    
    return os;
}

}
#endif
