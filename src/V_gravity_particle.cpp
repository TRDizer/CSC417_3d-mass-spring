#include <V_gravity_particle.h>

void V_gravity_particle(double &V, Eigen::Ref<const Eigen::Vector3d> q,  double mass, Eigen::Ref<const Eigen::Vector3d> g) {
    // V = mgh and assume g is directional gravitational acceleration (gx, gy, gz) which can potentially be (0, g, 0)
    // g dot q makes sense if the gravitational force is not downwards
    V = mass * std::abs(g.dot(q));
}