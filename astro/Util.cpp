#include "Util.h"

#include <cspice/SpiceUsr.h>

namespace astro {

double pi()
{
    return pi_c();
}

double degreesPerRadian()
{
    return dpr_c();
}

double radiansPerDegree()
{
    return 1.0 / dpr_c();
}

void raDecToVec(double range, double ra, double dec, Vec3* res)
{
    radrec_c(range, ra, dec, &(res->x));
}

void vecToRaDec(const Vec3& pos, double* range, double* ra, double* dec)
{
    recrad_c(const_cast<double*>(&(pos.x)), range, ra, dec);
}

double wrap(double value, double start, double end)
{
    const double width       = end - start;
    const double offsetValue = value - start;
    return (offsetValue - (std::floor(offsetValue / width) * width)) + start;
}

Vec3 transform(const Quat& q1, const Vec3& v, const Quat& q2)
{
    // Pure quaternion from v
    Quat q_v(0.0, v.x, v.y, v.z);
    Quat res = q1 * q_v * q2;
    return Vec3(res.x, res.y, res.z);
}


std::ostream& operator<<(std::ostream& os, const Vec3& v)
{
    os << "[ " << v.x << ", " << v.y << ", " << v.z << " ]";
    return os;
}

std::ostream& operator<<(std::ostream& os, const Quat& q)
{
    os << "[ " << q.x << "i + " << q.y << "j + " << q.z << "k + " << q.w << " ]";
    return os;
}

std::ostream& operator<<(std::ostream& os, const Mat3& m)
{
    // GLM is column-major: m[col][row]. Print row by row.
    for (int row = 0; row < 3; ++row)
        os << "[ " << m[0][row] << ", " << m[1][row] << ", " << m[2][row] << " ]\n";
    return os;
}

Quat getProgradeOrientation(const Vec3& pos, const Vec3& velocity)
{
    Vec3 x = glm::normalize(velocity);
    Vec3 y = -glm::normalize(pos);
    Vec3 z = glm::normalize(glm::cross(x, y));
    // Orthogonalize
    y = glm::normalize(glm::cross(z, x));

    // Build rotation matrix (columns = basis vectors) and convert to quaternion
    Mat3 rot(x, y, z); // col 0=x, col 1=y, col 2=z
    return glm::quat_cast(rot);
}

Quat getRetrogradeOrientation(const Vec3& pos, const Vec3& velocity)
{
    Vec3 x = -glm::normalize(velocity);
    Vec3 y = glm::normalize(pos);
    Vec3 z = glm::normalize(glm::cross(x, y));
    // Orthogonalize
    y = glm::normalize(glm::cross(z, x));

    Mat3 rot(x, y, z);
    return glm::quat_cast(rot);
}

} // namespace astro
