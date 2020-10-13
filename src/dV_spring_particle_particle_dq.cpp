#include <dV_spring_particle_particle_dq.h>

using namespace std;

void dV_spring_particle_particle_dq(Eigen::Ref<Eigen::Vector6d> f, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d>     q1, double l0, double stiffness) {
    f.setZero();
    // double dV_constant = stiffness * (1 - l0 / (q1-q0).norm());
    // // f = dV(q)/dq = dV(q)/d(q0_x, q0_y, q0_z, q1_x, q1_y, q1_z)
    // // hence the the last term of q0_? would yield a -1 to flip the difference term
    // // f << dV_constant * (q0 - q1)(0),
    // //      dV_constant * (q0 - q1)(1),
    // //      dV_constant * (q0 - q1)(2),
    // //      dV_constant * (q1 - q0)(0),
    // //      dV_constant * (q1 - q0)(1),
    // //      dV_constant * (q1 - q0)(2);
    // f << dV_constant * (q0 - q1), 
    //      dV_constant * (q1 - q0);
    double q0x = q0(0);
    double q0y = q0(1);
    double q0z = q0(2);
    double q1x = q1(0);
    double q1y = q1(1);
    double q1z = q1(2);
    f(0) = (q0x*2.0-q1x*2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z))*(-1.0/2.0);
    f(1) = (q0y*2.0-q1y*2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z))*(-1.0/2.0);
    f(2) = (q0z*2.0-q1z*2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z))*(-1.0/2.0);
    f(3) = ((q0x*2.0-q1x*2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))/2.0;
    f(4) = ((q0y*2.0-q1y*2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))/2.0;
    f(5) = ((q0z*2.0-q1z*2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))/2.0;
    f *= stiffness;
}