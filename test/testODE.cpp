
#include "../astro/ODE.h"
#include "../astro/State.h"
#include "../astro/Util.h"
#include <gtest/gtest.h>

#include <cmath>

using namespace astro;

class ODETest : public ::testing::Test {

protected:
    ODETest();

    virtual ~ODETest();

    // Code here will be called immediately after the constructor (right
    // before each test).
    virtual void SetUp();

    // Code here will be called immediately after each test (right
    // before the destructor).
    virtual void TearDown();

    double mu_earth;
    astro::ODE  ode0;
};



ODETest::ODETest()
  :  mu_earth(398600.0)
{
    astro::Attractor earth = {Vec3(0.0), mu_earth};
    ode0.addAttractor(earth);



}

ODETest::~ODETest()
{

}

void ODETest::SetUp()
{
}

void ODETest::TearDown()
{
}

// ---------------------------------------------------------------------------
// ODE::setForce / setBodyForce / setMass tests
// ---------------------------------------------------------------------------

// With no applied force, acceleration equals gravity only.
TEST_F(ODETest, NoForceGravityOnly)
{
    // Circular orbit radius 7000 km, velocity = sqrt(mu/r)
    double r = 7000.0;
    double v = std::sqrt(mu_earth / r);
    PosState s;
    s.r = Vec3(r, 0.0, 0.0);
    s.v = Vec3(0.0, v, 0.0);

    EphemerisTime et(0.0);
    PosState sdot = ode0.rates(et, s);

    // Centripetal acceleration magnitude = mu/r^2
    double expected_a = mu_earth / (r * r);
    double actual_a   = glm::length(sdot.v);
    EXPECT_NEAR(actual_a, expected_a, 1.0e-9);
}

// setForce() adds a constant inertial-frame acceleration (f/m).
TEST_F(ODETest, InertialForceAddsAcceleration)
{
    double r = 7000.0;
    PosState s;
    s.r = Vec3(r, 0.0, 0.0);
    s.v = Vec3(0.0, 0.0, 0.0);

    EphemerisTime et(0.0);

    double mass     = 1000.0;               // kg
    Vec3   force    = Vec3(0.0, 100.0, 0.0); // 100 N in Y (inertial)
    double expected_extra = 100.0 / mass;   // 0.1 m/s^2

    ODE ode_f;
    ode_f.addAttractor({Vec3(0.0), mu_earth});
    ode_f.setMass(mass);
    ode_f.setForce(force);

    PosState sdot_f = ode_f.rates(et, s);
    PosState sdot_0 = ode0.rates(et, s);   // no force baseline

    // Y-component of acceleration should differ by exactly force/mass
    EXPECT_NEAR(sdot_f.v.y - sdot_0.v.y, expected_extra, 1.0e-12);
    // X-component (gravity only) should be unchanged
    EXPECT_NEAR(sdot_f.v.x, sdot_0.v.x, 1.0e-12);
}

// setBodyForce() with identity attitude equals setForce() in inertial frame.
TEST_F(ODETest, BodyForceIdentityAttitude)
{
    double r = 7000.0;
    PosState s;
    s.r = Vec3(r, 0.0, 0.0);
    s.v = Vec3(0.0, 0.0, 0.0);

    EphemerisTime et(0.0);

    double mass  = 500.0;
    Vec3   force = Vec3(10.0, 20.0, 30.0);
    Quat   identity(1.0, 0.0, 0.0, 0.0);

    ODE ode_body;
    ode_body.addAttractor({Vec3(0.0), mu_earth});
    ode_body.setMass(mass);
    ode_body.setBodyForce(force, identity);

    ODE ode_inertial;
    ode_inertial.addAttractor({Vec3(0.0), mu_earth});
    ode_inertial.setMass(mass);
    ode_inertial.setForce(force);

    PosState sdot_b = ode_body.rates(et, s);
    PosState sdot_i = ode_inertial.rates(et, s);

    EXPECT_NEAR(sdot_b.v.x, sdot_i.v.x, 1.0e-12);
    EXPECT_NEAR(sdot_b.v.y, sdot_i.v.y, 1.0e-12);
    EXPECT_NEAR(sdot_b.v.z, sdot_i.v.z, 1.0e-12);
}

// setBodyForce() with a 90° rotation about Z: body X → inertial Y.
TEST_F(ODETest, BodyForceRotated90DegZ)
{
    PosState s;
    s.r = Vec3(7000.0, 0.0, 0.0);
    s.v = Vec3(0.0, 0.0, 0.0);
    EphemerisTime et(0.0);

    double mass = 100.0;
    // 90° rotation about Z: body X maps to inertial Y
    Quat att = glm::angleAxis(glm::radians(90.0), Vec3(0.0, 0.0, 1.0));

    ODE ode_body;
    ode_body.addAttractor({Vec3(0.0), mu_earth});
    ode_body.setMass(mass);
    ode_body.setBodyForce(Vec3(mass, 0.0, 0.0), att); // 1 N/kg = 1 m/s^2 in body X

    PosState sdot_0 = ode0.rates(et, s);
    PosState sdot_b = ode_body.rates(et, s);

    // Body X force → inertial Y acceleration (1 m/s^2)
    EXPECT_NEAR(sdot_b.v.y - sdot_0.v.y,  1.0, 1.0e-10);
    // Body X force should contribute nothing to inertial X
    EXPECT_NEAR(sdot_b.v.x - sdot_0.v.x,  0.0, 1.0e-10);
    EXPECT_NEAR(sdot_b.v.z - sdot_0.v.z,  0.0, 1.0e-10);
}

// setForce(zero) after a non-zero force resets to gravity-only.
TEST_F(ODETest, ClearForce)
{
    PosState s;
    s.r = Vec3(7000.0, 0.0, 0.0);
    s.v = Vec3(0.0, 0.0, 0.0);
    EphemerisTime et(0.0);

    ODE ode_f;
    ode_f.addAttractor({Vec3(0.0), mu_earth});
    ode_f.setMass(200.0);
    ode_f.setForce(Vec3(0.0, 1000.0, 0.0));
    ode_f.setForce(Vec3(0.0));   // clear

    PosState sdot_f = ode_f.rates(et, s);
    PosState sdot_0 = ode0.rates(et, s);

    EXPECT_NEAR(sdot_f.v.x, sdot_0.v.x, 1.0e-12);
    EXPECT_NEAR(sdot_f.v.y, sdot_0.v.y, 1.0e-12);
    EXPECT_NEAR(sdot_f.v.z, sdot_0.v.z, 1.0e-12);
}

// ---------------------------------------------------------------------------

TEST_F(ODETest, RotODESingularInertia)
{
    Mat3 I_singular = Mat3(0.0);
    RotODE rode;
    ASSERT_THROW(rode.setInertialMatrix(I_singular), astro::AstroException);

    Mat3 I_id = Mat3(1.0);
    ASSERT_NO_THROW(rode.setInertialMatrix(I_id));

    // Test inertia set in construction of rot ode
    ASSERT_THROW(RotODE rode2(I_singular), AstroException);
    ASSERT_NO_THROW(RotODE rode2(I_id));




}


