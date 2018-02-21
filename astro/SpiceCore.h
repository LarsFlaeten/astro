#ifndef _INCLUDE_SPICE_CORE_H_
#define _INCLUDE_SPICE_CORE_H_

#include <mutex>
#include <vector>
#include <string>

#include "Exceptions.h"

namespace astro
{

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

private:
    // Global mutex for spice common resource access
    std::mutex  m;  

    // a list of currently loaded Spice Kernels
    std::vector<std::string>    loadedKernels;



};


// Singleton access to SpiceCore
SpiceCore&  Spice();



}
#endif
