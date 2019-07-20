#include "State.h"

#include "Exceptions.h"
#include "SpiceCore.h"
#include <cspice/SpiceUsr.h>

namespace astro
{
    PosState    PosState::transform(const ReferenceFrame& fromFr, const ReferenceFrame toFr, const EphemerisTime& et)
    {
        if(fromFr.getType() == ReferenceFrameType::BodyFixedRotating && toFr.getType() == ReferenceFrameType::BodyFixedRotating) {
            throw AstroException("Either to-frame or from-frame must be inertial when transforming state");

        }

        if(( fromFr.isJ2000() || fromFr.getType() == ReferenceFrameType::BodyFixedNonRotating ) &&
           ( toFr.isJ2000() || toFr.getType() == ReferenceFrameType::BodyFixedNonRotating ))
            return *this; // No state transform needed

        // Ok, so either to r from is inertial, and the other is body fixed. Nevertheless,
        // get the bodyfixed one
        bool fromInertial;
        int body;
        if(toFr.isJ2000() || toFr.getType() == ReferenceFrameType::BodyFixedNonRotating ) {
            fromInertial = false;
            body = fromFr.getCenterId();
        } else {
            fromInertial = true;
            body = toFr.getCenterId();
        }

        double tispm[6][6];
 
        // Lock spice and retrieve the rotation from inertial to body-fixed:
        {
            std::lock_guard<std::mutex> lock(Spice().mutex());
            tisbod_c("J2000", body, et.getETValue(), tispm);          
           
        }
            
        if(!fromInertial) {
            double tispm_inv[6][6];
            invstm_c(tispm, tispm_inv);
            for(int i = 0; i < 6; ++i)
                for(int j = 0; j < 6; ++j)
                    tispm[i][j] = tispm_inv[i][j];
        } 
 
        // TODO: We ca also retrieve the angular velocities from the transformation here
        double state[6];
        double tr_state[6];
        state[0] = r.x;
        state[1] = r.y;
        state[2] = r.z;
        state[3] = v.x;
        state[4] = v.y;
        state[5] = v.z;
       
        // Do the tranformation 
        mxvg_c(tispm, state, 6, 6, tr_state);

        PosState R;
        R.r.x = tr_state[0];
        R.r.y = tr_state[1];
        R.r.z = tr_state[2];
        R.v.x = tr_state[3];
        R.v.y = tr_state[4];
        R.v.z = tr_state[5];

        return R;








    }


std::ostream& operator << (std::ostream& os, const astro::PosState& s)
{
    os << "r: (" << s.r.x <<", " << s.r.y << ", " << s.r.z << ") [km]\n";
    os << "v: (" << s.v.x <<", " << s.v.y << ", " << s.v.z << ") [km/s]\n";
    return os;
}   


PosState operator/(const PosState& p1, const PosState& p2)
{
    return PosState( p1.r/p2.r , p1.r/p2.r);
}


PosState abs(const PosState& p)
{
    return PosState( vec3d(std::abs(p.r.x) , std::abs(p.r.y) , std::abs(p.r.z) ),
					vec3d(std::abs(p.v.x) , std::abs(p.v.y) , std::abs(p.v.z)));

}



}
