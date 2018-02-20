#include <iostream>

#include "test.h"


#include <cspice/SpiceUsr.h>


namespace astro
{


    void astro_func() {
        std::cout << "In astro::astrofunc..." << std::endl;

    }

    void load_kernel(const std::string& filename)
    {
        furnsh_c(filename.c_str());
    }
}
