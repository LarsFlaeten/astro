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

enum AbborationCorrection {
    None,               // No correction, i.e geometric state
    LightTime,          // Light Time corrected state
    LightTimeStellar,   // Light time + stellar abboration (apparent state)
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

    // Loads the given file into spices kernel pool
    void    loadKernel(const std::string& filename);

    // returns the relative geometric state beween two celestial objects
    void    getRelativeGeometricState(int id, int ref_id, const EphemerisTime& et, astro::PosState& state);


private:
    // Global mutex for spice common resource access
    std::mutex  m;  

    // a list of currently loaded Spice Kernels
    std::vector<std::string>    loadedKernels;

    // Returns the Spice code for abboration
    void    getAbborationCode(AbborationCorrection ac, std::string& code);

};


// Singleton access to SpiceCore
SpiceCore&  Spice();



}
#endif
