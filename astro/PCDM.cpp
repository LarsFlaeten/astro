#include "PCDM.h"
#include "Util.h"

namespace astro {

PCDM::Result PCDM::doStep(const RotODE& rode, const RotState& rs, const EphemerisTime& et, const TimeDelta& dt)
{
    double DT = dt.value;

    // Rotation at time n
    Quat qn     = rs.q;
    Quat qn_inv = glm::inverse(qn);

    // Angular velocity at time n as a pure quaternion (scalar = 0)
    // GLM dquat(w, x, y, z), so scalar=0 means (0, wx, wy, wz)
    Quat q_wn = Quat(0.0, rs.w.x, rs.w.y, rs.w.z);

    // Transform angular velocities to body frame
    Quat q_wbn = qn_inv * q_wn * qn;
    Vec3 wbn   = Vec3(q_wbn.x, q_wbn.y, q_wbn.z);

    // Derivatives at time n
    RotState rsn_dot = rode.rates(et, rs);

    // (57) Angular velocities at n+1/4 and n+1/2
    Vec3 wndot  = rsn_dot.w;
    Vec3 wbndot = transform(qn_inv, wndot, qn);
    Vec3 wbn14  = wbn + 0.25 * wbndot * DT;
    Vec3 wbn12  = wbn + 0.50 * wbndot * DT;

    // (58) Angular velocity in global frame at n+1/4
    Vec3 wn14 = transform(qn, wbn14, qn_inv);

    // (59) Predicted q'_n+1/2
    double wn14_l = glm::length(wn14);
    double F      = wn14_l * DT * 0.25;
    Vec3   tmp    = Vec3(0.0);
    if (wn14_l >= 1.0E-8)
        tmp = (wn14 / wn14_l) * std::sin(F);
    else
        F = 0.0;
    // GLM dquat(w, x, y, z)
    Quat qprime_n12 = Quat(std::cos(F), tmp.x, tmp.y, tmp.z) * qn;

    // (60) Angular velocity in global frame at n+1/2
    Vec3 wn12 = transform(qprime_n12, wbn12, glm::inverse(qprime_n12));

    // Derivatives at n+1/2
    RotState rprime_n12(qprime_n12, wn12);
    TimeDelta dt12(0.5 * DT);
    RotState rsn12_dot = rode.rates(et + dt12, rprime_n12);

    // (61) q at n+1
    double wn12_l = glm::length(wn12);
    F   = wn12_l * DT * 0.5;
    tmp = Vec3(0.0);
    if (wn12_l >= 1.0E-8)
        tmp = (wn12 / wn12_l) * std::sin(F);
    else
        F = 0.0;
    Quat qn1     = Quat(std::cos(F), tmp.x, tmp.y, tmp.z) * qn;
    Quat qn1_inv = glm::inverse(qn1);

    // (62) Angular velocities at n+1
    Vec3 wn12dot  = rsn12_dot.w;
    Vec3 wbn12dot = transform(qn1_inv, wn12dot, qn1);
    Vec3 wbn1     = wbn + wbn12dot * DT;

    // (63) Transform to global frame
    Vec3 wn1 = transform(qn1, wbn1, qn1_inv);

    Result res;
    res.rs.q = qn1;
    res.rs.w = wn1;
    res.et   = et + dt;
    return res;
}

std::vector<PCDM::Result> PCDM::doSteps(
    const RotODE& rode, const RotState& rs,
    const EphemerisTime& et0, const EphemerisTime& et1,
    const TimeDelta& dt)
{
    std::vector<Result> res;
    res.push_back({ rs, et0 });

    TimeDelta dti = dt;
    while (res.back().et < et1)
    {
        res.push_back(doStep(rode, res.back().rs, res.back().et, dti));
        if (res.back().et + dti > et1)
            dti = et1 - res.back().et;
    }

    return res;
}

} // namespace astro
