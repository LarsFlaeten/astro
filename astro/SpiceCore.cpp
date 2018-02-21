#include "SpiceCore.h"
#include "Exceptions.h"

#include <cspice/SpiceUsr.h>

namespace astro
{

// Singleton acces to SpiceCore:
SpiceCore&  Spice()
{
    static SpiceCore spice;
    return spice;
}


SpiceCore::SpiceCore()
{
    // We set spice error acting here

    // We change spice to not abort program, but return on error:
    char error[] = "RETURN";
    erract_c( "SET", 0 , error);

    // We set spice to not report errors, since we handle them ourselves
    char prt[] = "NONE";
    errprt_c (  "SET", 0, prt );
}

SpiceCore::~SpiceCore()
{

}

void SpiceCore::checkError()
{
    if(failed_c())
    {
        // Probably need the mutex here..
        std::lock_guard<std::mutex> lock(m);
        
        char msg[1841]; // 1840 is max lenght of Spice long message
        int lenout = 1841;

        getmsg_c("SHORT", lenout, msg);
        std::string short_msg(msg);

        getmsg_c("LONG", lenout, msg);
        std::string long_msg(msg);

        getmsg_c("EXPLAIN", lenout, msg);
        std::string short_msg_expl(msg);

        // reset Spice errors, and throw exception
        reset_c();

        throw SpiceException(short_msg, long_msg, short_msg_expl);

    }

}

std::mutex& SpiceCore::mutex()
{
    return m;
}

void    SpiceCore::loadKernel(const std::string& filename)
{
    {
        std::lock_guard<std::mutex> lock(m);
        furnsh_c(filename.c_str());
    }
    checkError();
    
    loadedKernels.push_back(filename);

}


}
