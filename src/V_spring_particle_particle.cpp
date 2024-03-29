#include <V_spring_particle_particle.h>
#include <cmath>
#include <iostream>

//the potential energy of a spring with 3D end points q0 and qd and undeformed length l0
void V_spring_particle_particle(double &V, Eigen ::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness) {

    V = 1/2. * stiffness * std::pow((q1 - q0).norm() - l0, 2);
    
}