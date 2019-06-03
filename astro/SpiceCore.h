#ifndef _INCLUDE_SPICE_CORE_H_
#define _INCLUDE_SPICE_CORE_H_

#include <mutex>
#include <vector>
#include <string>

#include "Exceptions.h"
#include "Time.h"
#include "State.h"

namespace astro
{

enum AberrationCorrection {
    None,               // No correction, i.e geometric state
    LightTime,          // Light Time corrected state
    LightTimeStellar,   // Light time + stellar aberration (apparent state)
    CNLightTime,        // Converged Newtonian LT
    CNLightTimeStellar  // Converged Newtonian LT + Stellar Abboration
};

class SpiceCore
{
public:
    SpiceCore();
    virtual ~SpiceCore();

    // Provides a global mutex to ensure thread safe access to Spice
    std::mutex&  mutex();

    // Check wether Spice has the error flag set, and throw exception if so.
    // The function resets Spice's error handling mechanism before throwing,
    // so it is safe to continue if an exception is caught and handled.
    // This method must be called after any single or sequential spice function calls.
    // It is tedious to have to call this method explcicitly, but since
    // Spice has no way of setting an error callback, this is assumed the
    // best choice for now..
    void    checkError();

    // Sets whether spice should print errors to screen, in addition to throwing exceptions
    void    reportErrors(bool rep);

    // Loads the given file into spices kernel pool
    void    loadKernel(const std::string& filename);

    // returns the relative geometric state beween two celestial objects
    // int tgt_id: Target id
    // int obs_id: Observer id
    // state: State to be written
    void    getRelativeGeometricState(int tgt_id, int obs_id, const EphemerisTime& et, astro::PosState& state);

    // returns the relative geometric position beween two celestial objects,
    // optionally corrected for light time and stellar abboration
    // int tgt_id: Target id
    // int obs_id: Observer id
    // state: State to be written
    void    getRelativePosition(int tgt_id, int obs_id, const EphemerisTime& et, mork::vec3d& pos, AberrationCorrection abcorr = None);


    // Returns a requested constant/constant set from Spice
    // typically "GM", "RADII" etc
    void    getPlanetaryConstants(int id, const std::string& item, int num, double* vals);



private:
    // Global mutex for spice common resource access
    std::mutex  m;  

    // a list of currently loaded Spice Kernels
    std::vector<std::string>    loadedKernels;

    // Returns the Spice code for abboration
    void    getAberrationCode(AberrationCorrection ac, std::string& code);


    bool    report_errors;
};

// Debug function
// Prints the loaded kernels from Spice to ostream
std::ostream& operator << (std::ostream& os, const SpiceCore& s);
 

// Singleton access to SpiceCore
SpiceCore&  Spice();



}
#endif
