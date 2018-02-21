#include <iostream>

#include <astro/SpiceCore.h>


int main(int argc, char **argv)
{
    std::cout << "In main..." << std::endl;


    astro::Spice().loadKernel("../data/spice/lsk/naif0012.tls");

    // This will throw an exception:
    astro::Spice().loadKernel("NoSuchFile");


    return 0;
}
