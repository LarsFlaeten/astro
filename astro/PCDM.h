#ifndef _ASTRO_PCDM_H_
#define _ASTRO_PCDM_H_

#include <vector>
#include "Time.h"
#include "State.h"
#include "ODE.h"


namespace astro
{

// Implementation of the Predictor-Corrector Direct Multiplication algorithm
// by Zhao/van Wachem for numerical integration of rotation ODEs with quternions
// http://calliope.dem.uniud.it/COST/downloads/paper1_BVW.pdf
class PCDM
{
public:
    struct Result
    {
        RotState rs;
        EphemerisTime et;
    };


    static Result doStep(const RotODE& rode, const RotState& s, const EphemerisTime& et, const TimeDelta& dt);            

    static std::vector<Result> doSteps(const RotODE& rode, const RotState& s, const EphemerisTime& et0, const EphemerisTime& et1, const TimeDelta& dt);            


};




}


#endif
