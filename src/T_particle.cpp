#include <T_particle.h>
#include <iostream>

void T_particle(double &T, Eigen::Ref<const Eigen::VectorXd> qdot, double mass) {

    // std::cout << "my qdot:\n" << qdot << std::endl;
    T = 1/2. * mass * qdot.dot(qdot);
}