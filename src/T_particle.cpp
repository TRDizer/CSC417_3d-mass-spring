#include <T_particle.h>

void T_particle(double &T, Eigen::Ref<const Eigen::VectorXd> qdot, double mass) {

    T = 1/2 * mass * qdot.dot(qdot);
}