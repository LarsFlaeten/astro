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
    report_errors = false;
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

void SpiceCore::reportErrors(bool rep)
{
    report_errors = rep;
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

        if(report_errors)
        {
            std::cout << "Spice ERROR:\n";
            std::cout << short_msg << "\n";
            std::cout << short_msg_expl << "\n";
            std::cout << long_msg << std::endl;
        }

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

void    SpiceCore::getRelativeGeometricState(int tgt_id, int obs_id, const EphemerisTime& et, astro::PosState& state)
{
    double s[6];
    double lt;
    {
        std::lock_guard<std::mutex> lock(m);
        if(obs_id == 0)
            spkssb_c(tgt_id, et.getETValue(), "J2000", s);
        else
            spkgeo_c(tgt_id, et.getETValue(), "J2000", obs_id, s, &lt);
    }
    checkError();
   
    state.r.x = s[0];
    state.r.y = s[1];
    state.r.z = s[2];
    state.v.x = s[3];
    state.v.y = s[4];
    state.v.z = s[5];






}

void    SpiceCore::getRelativePosition(int tgt_id, int obs_id, const EphemerisTime& et, mork::vec3d& pos, AberrationCorrection abcorr)
{
    double p[3];
    double lt;
    std::string ac;
    getAberrationCode(abcorr, ac);
    
    {
        std::lock_guard<std::mutex> lock(m);
        spkezp_c(tgt_id, et.getETValue(), "J2000", ac.c_str(), obs_id, p, &lt);
    }
    checkError();

    pos.x = p[0];
    pos.y = p[1];
    pos.z = p[2];
}

void    SpiceCore::getAberrationCode(AberrationCorrection ac, std::string& code)
{
    switch(ac) {
        case(None):
            code = "NONE";
            break;
        case(LightTime):
            code = "LT";
            break;
        case(LightTimeStellar):
            code = "LT+S";
            break;
        case(CNLightTime):
            code = "CN";
            break;
        case(CNLightTimeStellar):
            code = "CN+S";
            break;
        default:
            throw SpiceException("Abboration type not implemented", "", "");
    }
}

void   SpiceCore::getPlanetaryConstants(int id, const std::string& item, int num, double* vals)
{
    {
        std::lock_guard<std::mutex> lock(m);
        int read_num;
        auto s = std::to_string(id);
        bodvrd_c(s.c_str(), item.c_str(), num, &read_num, vals);
    }
    checkError();

}

// Prints various information about spice core
std::ostream& operator << (std::ostream& os, const SpiceCore& s)
{
    // Return the current number of kernels that have been loaded 
    // via the KEEPER interface that are of a specified type.  
    int count;
    ktotal_c("ALL", &count);
    
    os << "Kernel summary\n";
    os << "============================\n";

    if(count == 0)
        os << "No kernels loaded\n";
    else
        os << count << " kernel files loaded:\n";

    const int FILLEN =   256;
    const int TYPLEN =   33;
    const int SRCLEN =  256;

    int        handle;

    char      file  [FILLEN];
    char       filtyp[TYPLEN];
    char       source[SRCLEN];
    int    found;

    for(int i = 0; i < count; ++i)
    {
       kdata_c ( i,  "all",    FILLEN,   TYPLEN, SRCLEN, 
                                   file,   filtyp,  source,  &handle,  &found );  

       kinfo_c ( file, TYPLEN, SRCLEN, filtyp, source, &handle, &found );

        os << "  " << i+1 << ": " << file << "[" << filtyp << "], " << source << "\n";
    }

    return os;
}
 

}
