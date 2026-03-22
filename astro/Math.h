#pragma once

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtx/norm.hpp>

namespace astro {
    using Vec3 = glm::dvec3;
    using Vec4 = glm::dvec4;
    using Quat = glm::dquat;
    using Mat3 = glm::dmat3;
    using Mat4 = glm::dmat4;
} // namespace astro
