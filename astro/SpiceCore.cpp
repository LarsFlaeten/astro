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

void    SpiceCore::getRelativeGeometricState(int tgt_id, int obs_id, const EphemerisTime& et, astro::PosState& state, const ReferenceFrame& rf)
{
    double s[6];
    double lt;
    {
        std::lock_guard<std::mutex> lock(m);
        if(obs_id == 0)
            spkssb_c(tgt_id, et.getETValue(), rf.getName().c_str(), s);
        else
            spkgeo_c(tgt_id, et.getETValue(), rf.getName().c_str(), obs_id, s, &lt);
    }
    checkError();
   
    state.r.x = s[0];
    state.r.y = s[1];
    state.r.z = s[2];
    state.v.x = s[3];
    state.v.y = s[4];
    state.v.z = s[5];






}
    
void    SpiceCore::getRelativeState(int tgt_id, const Observer& obs, const EphemerisTime& et, astro::PosState& state, AberrationCorrection abcorr) {
    double s[6];
    double lt;
    std::string tgt_id_s = std::to_string(tgt_id);
    std::string obs_ref_fr_s = obs.getReferenceFrame().getName();
    std::string ac;
    getAberrationCode(abcorr, ac);
    {
        std::lock_guard<std::mutex> lock(m);
 
        // Get the observer state:
        PosState obsState = obs.getState();
        int centerObj = obs.getCenterObject();
        // Corrrect the state to SSB relative if the RF is J2000 (Since that is how spice works
        // when using J2000 in geometry calcs (SSB is the center))
        // In astro, we use the frame for rotation and allow J2000 to be used for other centerObjects
        // TODO: Clean up in all this, this might get messy....
        PosState c_state;
        if(obs.getReferenceFrame().isJ2000() && centerObj != 0) 
        {
            spkssb_c(centerObj, et.getETValue(), "J2000", (double*)&c_state);
            obsState += c_state; 
            centerObj = 0;
        }

        std::string c_obj_s = std::to_string(centerObj);
        // We can now use Spices observer state function:
       
        // Call spice and assign result to state
        spkcvo_c(tgt_id_s.c_str(), et.getETValue(), obs_ref_fr_s.c_str(), "OBSERVER",
            ac.c_str(), &obsState, et.getETValue(), c_obj_s.c_str(), obs_ref_fr_s.c_str(), s, &lt);
        
        state.r.x = s[0];
        state.r.y = s[1];
        state.r.z = s[2];
        state.v.x = s[3];
        state.v.y = s[4];
        state.v.z = s[5];


    }
    checkError();
   


}


void    SpiceCore::getRelativePosition(int tgt_id, int obs_id, const EphemerisTime& et, mork::vec3d& pos, AberrationCorrection abcorr, const ReferenceFrame& rf)
{
    double p[3];
    double lt;
    std::string ac;
    getAberrationCode(abcorr, ac);
    
    {
        std::lock_guard<std::mutex> lock(m);
        spkezp_c(tgt_id, et.getETValue(), rf.getName().c_str(), ac.c_str(), obs_id, p, &lt);
    }
    checkError();

    pos.x = p[0];
    pos.y = p[1];
    pos.z = p[2];
}

void    SpiceCore::getRelativePosition(int tgt_id, const Observer& obs, const EphemerisTime& et, mork::vec3d& pos, AberrationCorrection abcorr)
{
    PosState tgt_state;
    getRelativeState(tgt_id, obs, et, tgt_state, abcorr);
    pos = tgt_state.r;
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

void   SpiceCore::getPlanetaryConstants(int id, const std::string& item, double& val)
{
    getPlanetaryConstants(id, item, 1, &val);
}

void   SpiceCore::getPlanetaryConstants(int id, const std::string& item, mork::vec3d& val)
{
    getPlanetaryConstants(id, item, 3, &(val.x));
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
