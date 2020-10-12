#include <assemble_stiffness.h>

void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double k) { 
    //     d(q0x)(q0x)    d(q0x)(q0y)    d(q0x)(q0z)    d(q0x)(q1x)    d(q0x)(q1y)    d(q0x)(q1z)    
    //     d(q0y)(q0x)    d(q0y)(q0y)    d(q0y)(q0z)    d(q0y)(q1x)    d(q0y)(q1y)    d(q0y)(q1z)    
    //     d(q0z)(q0x)    d(q0z)(q0y)    d(q0z)(q0z)    d(q0z)(q1x)    d(q0z)(q1y)    d(q0z)(q1z)    
    // H = 
    //     d(q1x)(q0x)    d(q1x)(q0y)    d(q1x)(q0z)    d(q1x)(q1x)    d(q1x)(q1y)    d(q1x)(q1z)    
    //     d(q1y)(q0x)    d(q1y)(q0y)    d(q1y)(q0z)    d(q1y)(q1x)    d(q1y)(q1y)    d(q1y)(q1z)    
    //     d(q1z)(q0x)    d(q1z)(q0y)    d(q1z)(q0z)    d(q1z)(q1x)    d(q1z)(q1y)    d(q1z)(q1z)    
    //
    //     H(q0/q0)    H(q0/q1)    ........    H(q0/qn)
    //     H(q1/q0)    ........    ........    ........
    //     ........    ........    ........    ........
    // K = ........    ........    ........    ........
    //     ........    ........    ........    ........
    //     ........    ........    ........    ........
    //     H(qn/q0)    ........    ........    H(qn/qn)
    K.resize(q.size(), q.size());
    K.setZero();
    Eigen::Matrix66d spring_i_H;
    std::vector<Eigen::Triplet<double>> K_entries;
    int q0_index, q1_index;
    for (int spring_i = 0; spring_i < E.rows(); spring_i++) {
        q0_index = 3 * E(spring_i, 0);
        q1_index = 3 * E(spring_i, 1);
        Eigen::Vector3d q0 = q.segment<3>(q0_index);
        Eigen::Vector3d q1 = q.segment<3>(q1_index);
        d2V_spring_particle_particle_dq2(spring_i_H, q0, q1, l0(spring_i), k);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                K_entries.push_back(Eigen::Triplet<double>(q0_index + i, q0_index + j, -spring_i_H(i,j)));
                K_entries.push_back(Eigen::Triplet<double>(q0_index + i, q1_index + j, -spring_i_H(i,j + 3)));
                K_entries.push_back(Eigen::Triplet<double>(q1_index + i, q0_index + j, -spring_i_H(i + 3,j)));
                K_entries.push_back(Eigen::Triplet<double>(q1_index + i, q1_index + j, -spring_i_H(i + 3,j + 3)));
            }
        }
    }
    K.setFromTriplets(K_entries.begin(), K_entries.end());
};