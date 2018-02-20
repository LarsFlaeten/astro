#include <iostream>

#include <astro/test.h>


int main()
{
    std::cout << "In main..." << std::endl;

    astro::astro_func();

    astro::load_kernel("naif0012.tls");


    return 0;
}
