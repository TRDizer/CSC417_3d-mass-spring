#include <dV_gravity_particle_dq.h>

void dV_gravity_particle_dq(Eigen::Ref<Eigen::Vector3d> f,  double mass, Eigen::Ref<const Eigen::Vector3d> g) {
    // Gravitational force should still be directional despite the energy is abs'ed
    f = mass * g;
}