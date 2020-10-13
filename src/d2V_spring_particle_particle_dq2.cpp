#include <d2V_spring_particle_particle_dq2.h>
#include <cmath>
#include <iostream>

#define GET_dq_SIGN(i,k)                     (i == 1 ? 1 : -1)
#define GET_dq_DIFF_TERM(diff_term, i, k)    (diff_term(k))

using namespace std;

void d2V_spring_particle_particle_dq2(Eigen::Ref<Eigen::Matrix66d> H, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness) {
    // dV = stiffness * dV_term0 * (GET_dq_SIGN(i,k) where ik in dq_ik) * (GET_dq_DIFF_TERM(diff_term, i, k) where diff_term = q1-q0, ik in dq_ik) 
    // const Eigen::Vector3d dx = q1 - q0;
    // double dV_term0 = 1 - l0 / dx.norm();
    // auto dV_term1 = [=] (int i, int k) {
    //     return GET_dq_SIGN(i,k) * GET_dq_DIFF_TERM(dx, i, k);
    // };
    // //     d(q0x)(q0x)    d(q0x)(q0y)    d(q0x)(q0z)    d(q0x)(q1x)    d(q0x)(q1y)    d(q0x)(q1z)    
    // //     d(q0y)(q0x)    d(q0y)(q0y)    d(q0y)(q0z)    d(q0y)(q1x)    d(q0y)(q1y)    d(q0y)(q1z)    
    // //     d(q0z)(q0x)    d(q0z)(q0y)    d(q0z)(q0z)    d(q0z)(q1x)    d(q0z)(q1y)    d(q0z)(q1z)    
    // // H = 
    // //     d(q1x)(q0x)    d(q1x)(q0y)    d(q1x)(q0z)    d(q1x)(q1x)    d(q1x)(q1y)    d(q1x)(q1z)    
    // //     d(q1y)(q0x)    d(q1y)(q0y)    d(q1y)(q0z)    d(q1y)(q1x)    d(q1y)(q1y)    d(q1y)(q1z)    
    // //     d(q1z)(q0x)    d(q1z)(q0y)    d(q1z)(q0z)    d(q1z)(q1x)    d(q1z)(q1y)    d(q1z)(q1z)    
    // // product rule, here we go (tired face)
    // auto dV_term0_dqik = [=] (int i, int k) {
    //     return l0 / std::pow(dx.norm(), 3) * GET_dq_SIGN(i,k) * GET_dq_DIFF_TERM(dx, i, k);
    // };
    // auto dV_term1nm_dqik = [=] (int n, int m, int i, int k) {
    //     return (m == k) ? GET_dq_SIGN(n,m) * GET_dq_SIGN(i,k): 0;
    // };
    // // d(q_ab)(q_cd)
    // int a,b,c,d;
    // for(int H_j = 0; H_j < 6; H_j++) {
    //     a = (int) (H_j / 3);
    //     b = H_j % 3;
    //     for (int H_k = 0; H_k <6; H_k++) {
    //         c = (int) (H_k / 3);
    //         d = H_k % 3;
    //         H(H_j, H_k) = stiffness * (
    //             dV_term0_dqik(c, d) * dV_term1(a,b) +
    //             dV_term0 * dV_term1nm_dqik(a,b,c,d)
    //         );
    //     }
    // }
    // ===================================================================================================
    // d2V/dq2 = k * (l0 / pow(l,3) * BMq * BMq) + BM * k (1 - l0 / l)
    //        I    -I
    //  BM = 
    //       -I     I
    // double term1 = l0 / dx.dot(dx) / dx.norm();
    // double term2 = 1 - l0 / dx.norm();
    // Eigen::MatrixXd BMqBMq(6,6);
    // BMqBMq << -dx * -dx.transpose(), -dx * dx.transpose(),
    //            dx * -dx.transpose(), dx * dx.transpose();
    // for(int H_j = 0; H_j < 3; H_j++) {
    //     for (int H_k = 0; H_k <3; H_k++) {
    //         H(0 + H_j, 0 + H_k) = term1 * BMqBMq(0 + H_j, 0 + H_k) + (H_j == H_k ? term2 : 0);
    //         H(0 + H_j, 3 + H_k) = term1 * BMqBMq(0 + H_j, 0 + H_k) + (H_j == H_k ? -term2 : 0);
    //         H(3 + H_j, 0 + H_k) = term1 * BMqBMq(0 + H_j, 0 + H_k) + (H_j == H_k ? -term2 : 0);
    //         H(3 + H_j, 3 + H_k) = term1 * BMqBMq(0 + H_j, 0 + H_k) + (H_j == H_k ? term2 : 0);
    //     }
    // }
    // H *= stiffness;
    double q0x = q0(0);
    double q0y = q0(1);
    double q0z = q0(2);
    double q1x = q1(0);
    double q1y = q1(1);
    double q1z = q1(2);
    H(0,0) = pow(q0x*2.0-q1x*2.0,2.0)/(q0x*(q0x-q1x)*4.0+q0y*(q0y-q1y)*4.0-q1x*(q0x-q1x)*4.0+q0z*(q0z-q1z)*4.0-q1y*(q0y-q1y)*4.0-q1z*(q0z-q1z)*4.0)-(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z))+(pow(q0x*2.0-q1x*2.0,2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(0,1) = ((q0x*2.0-q1x*2.0)*(q0y*2.0-q1y*2.0))/(q0x*(q0x-q1x)*4.0+q0y*(q0y-q1y)*4.0-q1x*(q0x-q1x)*4.0+q0z*(q0z-q1z)*4.0-q1y*(q0y-q1y)*4.0-q1z*(q0z-q1z)*4.0)+((q0x*2.0-q1x*2.0)*(q0y*2.0-q1y*2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(0,2) = ((q0x*2.0-q1x*2.0)*(q0z*2.0-q1z*2.0))/(q0x*(q0x-q1x)*4.0+q0y*(q0y-q1y)*4.0-q1x*(q0x-q1x)*4.0+q0z*(q0z-q1z)*4.0-q1y*(q0y-q1y)*4.0-q1z*(q0z-q1z)*4.0)+((q0x*2.0-q1x*2.0)*(q0z*2.0-q1z*2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(0,3) = (pow(q0x*2.0-q1x*2.0,2.0)*(-1.0/4.0))/(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z))+(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z))-(pow(q0x*2.0-q1x*2.0,2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(0,4) = ((q0x*2.0-q1x*2.0)*(q0y*2.0-q1y*2.0)*(-1.0/4.0))/(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z))-((q0x*2.0-q1x*2.0)*(q0y*2.0-q1y*2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(0,5) = ((q0x*2.0-q1x*2.0)*(q0z*2.0-q1z*2.0)*(-1.0/4.0))/(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z))-((q0x*2.0-q1x*2.0)*(q0z*2.0-q1z*2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(1,0) = ((q0x*2.0-q1x*2.0)*(q0y*2.0-q1y*2.0))/(q0x*(q0x-q1x)*4.0+q0y*(q0y-q1y)*4.0-q1x*(q0x-q1x)*4.0+q0z*(q0z-q1z)*4.0-q1y*(q0y-q1y)*4.0-q1z*(q0z-q1z)*4.0)+((q0x*2.0-q1x*2.0)*(q0y*2.0-q1y*2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(1,1) = pow(q0y*2.0-q1y*2.0,2.0)/(q0x*(q0x-q1x)*4.0+q0y*(q0y-q1y)*4.0-q1x*(q0x-q1x)*4.0+q0z*(q0z-q1z)*4.0-q1y*(q0y-q1y)*4.0-q1z*(q0z-q1z)*4.0)-(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z))+(pow(q0y*2.0-q1y*2.0,2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(1,2) = ((q0y*2.0-q1y*2.0)*(q0z*2.0-q1z*2.0))/(q0x*(q0x-q1x)*4.0+q0y*(q0y-q1y)*4.0-q1x*(q0x-q1x)*4.0+q0z*(q0z-q1z)*4.0-q1y*(q0y-q1y)*4.0-q1z*(q0z-q1z)*4.0)+((q0y*2.0-q1y*2.0)*(q0z*2.0-q1z*2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(1,3) = ((q0x*2.0-q1x*2.0)*(q0y*2.0-q1y*2.0)*(-1.0/4.0))/(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z))-((q0x*2.0-q1x*2.0)*(q0y*2.0-q1y*2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(1,4) = (pow(q0y*2.0-q1y*2.0,2.0)*(-1.0/4.0))/(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z))+(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z))-(pow(q0y*2.0-q1y*2.0,2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(1,5) = ((q0y*2.0-q1y*2.0)*(q0z*2.0-q1z*2.0)*(-1.0/4.0))/(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z))-((q0y*2.0-q1y*2.0)*(q0z*2.0-q1z*2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(2,0) = ((q0x*2.0-q1x*2.0)*(q0z*2.0-q1z*2.0))/(q0x*(q0x-q1x)*4.0+q0y*(q0y-q1y)*4.0-q1x*(q0x-q1x)*4.0+q0z*(q0z-q1z)*4.0-q1y*(q0y-q1y)*4.0-q1z*(q0z-q1z)*4.0)+((q0x*2.0-q1x*2.0)*(q0z*2.0-q1z*2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(2,1) = ((q0y*2.0-q1y*2.0)*(q0z*2.0-q1z*2.0))/(q0x*(q0x-q1x)*4.0+q0y*(q0y-q1y)*4.0-q1x*(q0x-q1x)*4.0+q0z*(q0z-q1z)*4.0-q1y*(q0y-q1y)*4.0-q1z*(q0z-q1z)*4.0)+((q0y*2.0-q1y*2.0)*(q0z*2.0-q1z*2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(2,2) = pow(q0z*2.0-q1z*2.0,2.0)/(q0x*(q0x-q1x)*4.0+q0y*(q0y-q1y)*4.0-q1x*(q0x-q1x)*4.0+q0z*(q0z-q1z)*4.0-q1y*(q0y-q1y)*4.0-q1z*(q0z-q1z)*4.0)-(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z))+(pow(q0z*2.0-q1z*2.0,2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(2,3) = ((q0x*2.0-q1x*2.0)*(q0z*2.0-q1z*2.0)*(-1.0/4.0))/(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z))-((q0x*2.0-q1x*2.0)*(q0z*2.0-q1z*2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(2,4) = ((q0y*2.0-q1y*2.0)*(q0z*2.0-q1z*2.0)*(-1.0/4.0))/(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z))-((q0y*2.0-q1y*2.0)*(q0z*2.0-q1z*2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(2,5) = (pow(q0z*2.0-q1z*2.0,2.0)*(-1.0/4.0))/(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z))+(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z))-(pow(q0z*2.0-q1z*2.0,2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(3,0) = (pow(q0x*2.0-q1x*2.0,2.0)*(-1.0/4.0))/(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z))+(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z))-(pow(q0x*2.0-q1x*2.0,2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(3,1) = ((q0x*2.0-q1x*2.0)*(q0y*2.0-q1y*2.0)*(-1.0/4.0))/(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z))-((q0x*2.0-q1x*2.0)*(q0y*2.0-q1y*2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(3,2) = ((q0x*2.0-q1x*2.0)*(q0z*2.0-q1z*2.0)*(-1.0/4.0))/(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z))-((q0x*2.0-q1x*2.0)*(q0z*2.0-q1z*2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(3,3) = pow(q0x*2.0-q1x*2.0,2.0)/(q0x*(q0x-q1x)*4.0+q0y*(q0y-q1y)*4.0-q1x*(q0x-q1x)*4.0+q0z*(q0z-q1z)*4.0-q1y*(q0y-q1y)*4.0-q1z*(q0z-q1z)*4.0)-(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z))+(pow(q0x*2.0-q1x*2.0,2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(3,4) = ((q0x*2.0-q1x*2.0)*(q0y*2.0-q1y*2.0))/(q0x*(q0x-q1x)*4.0+q0y*(q0y-q1y)*4.0-q1x*(q0x-q1x)*4.0+q0z*(q0z-q1z)*4.0-q1y*(q0y-q1y)*4.0-q1z*(q0z-q1z)*4.0)+((q0x*2.0-q1x*2.0)*(q0y*2.0-q1y*2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(3,5) = ((q0x*2.0-q1x*2.0)*(q0z*2.0-q1z*2.0))/(q0x*(q0x-q1x)*4.0+q0y*(q0y-q1y)*4.0-q1x*(q0x-q1x)*4.0+q0z*(q0z-q1z)*4.0-q1y*(q0y-q1y)*4.0-q1z*(q0z-q1z)*4.0)+((q0x*2.0-q1x*2.0)*(q0z*2.0-q1z*2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(4,0) = ((q0x*2.0-q1x*2.0)*(q0y*2.0-q1y*2.0)*(-1.0/4.0))/(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z))-((q0x*2.0-q1x*2.0)*(q0y*2.0-q1y*2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(4,1) = (pow(q0y*2.0-q1y*2.0,2.0)*(-1.0/4.0))/(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z))+(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z))-(pow(q0y*2.0-q1y*2.0,2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(4,2) = ((q0y*2.0-q1y*2.0)*(q0z*2.0-q1z*2.0)*(-1.0/4.0))/(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z))-((q0y*2.0-q1y*2.0)*(q0z*2.0-q1z*2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(4,3) = ((q0x*2.0-q1x*2.0)*(q0y*2.0-q1y*2.0))/(q0x*(q0x-q1x)*4.0+q0y*(q0y-q1y)*4.0-q1x*(q0x-q1x)*4.0+q0z*(q0z-q1z)*4.0-q1y*(q0y-q1y)*4.0-q1z*(q0z-q1z)*4.0)+((q0x*2.0-q1x*2.0)*(q0y*2.0-q1y*2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(4,4) = pow(q0y*2.0-q1y*2.0,2.0)/(q0x*(q0x-q1x)*4.0+q0y*(q0y-q1y)*4.0-q1x*(q0x-q1x)*4.0+q0z*(q0z-q1z)*4.0-q1y*(q0y-q1y)*4.0-q1z*(q0z-q1z)*4.0)-(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z))+(pow(q0y*2.0-q1y*2.0,2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(4,5) = ((q0y*2.0-q1y*2.0)*(q0z*2.0-q1z*2.0))/(q0x*(q0x-q1x)*4.0+q0y*(q0y-q1y)*4.0-q1x*(q0x-q1x)*4.0+q0z*(q0z-q1z)*4.0-q1y*(q0y-q1y)*4.0-q1z*(q0z-q1z)*4.0)+((q0y*2.0-q1y*2.0)*(q0z*2.0-q1z*2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(5,0) = ((q0x*2.0-q1x*2.0)*(q0z*2.0-q1z*2.0)*(-1.0/4.0))/(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z))-((q0x*2.0-q1x*2.0)*(q0z*2.0-q1z*2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(5,1) = ((q0y*2.0-q1y*2.0)*(q0z*2.0-q1z*2.0)*(-1.0/4.0))/(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z))-((q0y*2.0-q1y*2.0)*(q0z*2.0-q1z*2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(5,2) = (pow(q0z*2.0-q1z*2.0,2.0)*(-1.0/4.0))/(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z))+(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z))-(pow(q0z*2.0-q1z*2.0,2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(5,3) = ((q0x*2.0-q1x*2.0)*(q0z*2.0-q1z*2.0))/(q0x*(q0x-q1x)*4.0+q0y*(q0y-q1y)*4.0-q1x*(q0x-q1x)*4.0+q0z*(q0z-q1z)*4.0-q1y*(q0y-q1y)*4.0-q1z*(q0z-q1z)*4.0)+((q0x*2.0-q1x*2.0)*(q0z*2.0-q1z*2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(5,4) = ((q0y*2.0-q1y*2.0)*(q0z*2.0-q1z*2.0))/(q0x*(q0x-q1x)*4.0+q0y*(q0y-q1y)*4.0-q1x*(q0x-q1x)*4.0+q0z*(q0z-q1z)*4.0-q1y*(q0y-q1y)*4.0-q1z*(q0z-q1z)*4.0)+((q0y*2.0-q1y*2.0)*(q0z*2.0-q1z*2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H(5,5) = pow(q0z*2.0-q1z*2.0,2.0)/(q0x*(q0x-q1x)*4.0+q0y*(q0y-q1y)*4.0-q1x*(q0x-q1x)*4.0+q0z*(q0z-q1z)*4.0-q1y*(q0y-q1y)*4.0-q1z*(q0z-q1z)*4.0)-(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z))+(pow(q0z*2.0-q1z*2.0,2.0)*(l0-sqrt(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z)))*1.0/pow(q0x*(q0x-q1x)+q0y*(q0y-q1y)-q1x*(q0x-q1x)+q0z*(q0z-q1z)-q1y*(q0y-q1y)-q1z*(q0z-q1z),3.0/2.0))/4.0;
    H *= stiffness;
}