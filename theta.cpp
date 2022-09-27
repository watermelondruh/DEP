#include "theta.h"

Theta::Theta(Vector3d x0, Vector3d x1, Vector3d x2, Vector3d x3) {
    m_e0 = x1 - x0;
    m_e1 = x2 - x0;
    m_e2 = x3 - x0;
    m_e3 = x2 - x1;
    m_e4 = x3 - x1;

    m_cosA1 = m_e0.dot(m_e1) / (m_e0.norm() * m_e1.norm());
    m_cosA2 = m_e0.dot(m_e2) / (m_e0.norm() * m_e2.norm());
    m_cosA3 = m_e0.dot(m_e3) / (m_e0.norm() * m_e3.norm());
    m_cosA4 = m_e0.dot(m_e4) / (m_e0.norm() * m_e4.norm());

    m_sinA1 = (m_e0.cross(m_e1) / (m_e0.norm() * m_e1.norm())).norm();
    m_sinA2 = (m_e0.cross(m_e2) / (m_e0.norm() * m_e2.norm())).norm();
    m_sinA3 = - (m_e0.cross(m_e3) / (m_e0.norm() * m_e3.norm())).norm();
    m_sinA4 = - (m_e0.cross(m_e4) / (m_e0.norm() * m_e4.norm())).norm();

    m_nn1 = m_e0.cross(m_e3);
    m_nn1 = m_nn1 / m_nn1.norm();
    m_nn2 = m_e0.cross(m_e4);
    m_nn2 = m_nn2 / m_nn2.norm();

    m_h1 = m_e0.norm() * m_sinA1;
    m_h2 = m_e0.norm() * m_sinA2;
    m_h3 = - m_e0.norm() * m_sinA3;
    m_h4 = - m_e0.norm() * m_sinA4;
    m_h01 = m_e1.norm() * m_sinA1;
    m_h02 = m_e2.norm() * m_sinA2;
}

Theta::~Theta() { 
    ;
}

double Theta::getTheta() {
    Vector3d n0 = m_e0.cross(m_e1);
    Vector3d n1 = m_e2.cross(m_e0);

    Vector3d w = n0.cross(n1);
    double theta = atan(w.norm() / (n0.dot(n1)));
    if (m_e0.dot(w) < 0)
        theta = -theta;
    
    return theta;
}

VectorXd Theta::gradTheta() {
    VectorXd grad = VectorXd::Zero(12);
    grad.segment(0, 3) = m_cosA3 * m_nn1 / m_h3 + m_cosA4 * m_nn2 / m_h4;
    grad.segment(3, 3) = m_cosA1 * m_nn1 / m_h1 + m_cosA2 * m_nn2 / m_h2;
    grad.segment(6, 3) = - m_nn1 / m_h01;
    grad.segment(9, 3) = - m_nn2 / m_h02;

    return grad;
}

MatrixXd Theta::hessTheta() {
    VectorXd grad = gradTheta();

    Vector3d m_m1 = m_nn1.cross(m_e1) / m_e1.norm();
    Vector3d m_m2 = - m_nn2.cross(m_e2) / m_e2.norm();
    Vector3d m_m3 = - m_nn1.cross(m_e3) / m_e3.norm();
    Vector3d m_m4 = m_nn2.cross(m_e4) / m_e4.norm();
    Vector3d m_m01 = - m_nn1.cross(m_e0) / m_e0.norm();
    Vector3d m_m02 = m_nn2.cross(m_e0) / m_e0.norm();

    Matrix3d M331 = m_cosA3 / (m_h3 * m_h3) * m_m3 * m_nn1.transpose();
    Matrix3d M311 = m_cosA3 / (m_h3 * m_h1) * m_m1 * m_nn1.transpose();
    Matrix3d M131 = m_cosA1 / (m_h1 * m_h3) * m_m3 * m_nn1.transpose();
    Matrix3d M3011 = m_cosA3 / (m_h3 * m_h01) * m_m01 * m_nn1.transpose();
    Matrix3d M111 = m_cosA1 / (m_h1 * m_h1) * m_m1 * m_nn1.transpose();
    Matrix3d M1011 = m_cosA1 / (m_h1 * m_h01) * m_m01 * m_nn1.transpose();

    Matrix3d M442 = m_cosA4 / (m_h4 * m_h4) * m_m4 * m_nn2.transpose();
    Matrix3d M422 = m_cosA4 / (m_h4 * m_h2) * m_m2 * m_nn2.transpose();
    Matrix3d M242 = m_cosA2 / (m_h2 * m_h4) * m_m4 * m_nn2.transpose();
    Matrix3d M4022 = m_cosA4 / (m_h4 * m_h02) * m_m02 * m_nn2.transpose();
    Matrix3d M222 = m_cosA2 / (m_h2 * m_h2) * m_m2 * m_nn2.transpose();
    Matrix3d M2022 = m_cosA2 / (m_h2 * m_h02) * m_m02 * m_nn2.transpose();

    Matrix3d B1 = 1 / (m_e0.dot(m_e0)) * m_nn1 * m_m01.transpose();
    Matrix3d B2 = 1 / (m_e0.dot(m_e0)) * m_nn2 * m_m02.transpose();

    Matrix3d N13 = 1 / (m_h01 * m_h3) * m_nn1 * m_m3.transpose();
    Matrix3d N24 = 1 / (m_h02 * m_h4) * m_nn2 * m_m4.transpose();
    Matrix3d N11 = 1 / (m_h01 * m_h1) * m_nn1 * m_m1.transpose();
    Matrix3d N22 = 1 / (m_h02 * m_h2) * m_nn2 * m_m2.transpose();
    Matrix3d N101 = 1 / (m_h01 * m_h01) * m_nn1 * m_m01.transpose();
    Matrix3d N202 = 1 / (m_h02 * m_h02) * m_nn2 * m_m02.transpose();

    MatrixXd hess = MatrixXd::Zero(12, 12);
    hess.block(0, 0, 3, 3) = M331 + M331.transpose() - B1 + M442 + M442.transpose() - B2;
    hess.block(0, 3, 3, 3) = M331 + M131.transpose() + B1 + M422 + M242.transpose() + B2;
    hess.block(0, 6, 3, 3) = M3011 - N13;
    hess.block(0, 9, 3, 3) = M4022 - N24;
    hess.block(3, 3, 3, 3) = M111 + M111.transpose() - B1 + M222 + M222.transpose() - B2;
    hess.block(3, 6, 3, 3) = M1011 - N11;
    hess.block(3, 9, 3, 3) = M2022 - N22;
    hess.block(6, 6, 3, 3) = - N101 - N101.transpose();
    hess.block(9, 9, 3, 3) = - N202 - N202.transpose();
    hess.block(3, 0, 3, 3) = hess.block(0, 3, 3, 3).transpose();
    hess.block(6, 0, 3, 3) = hess.block(0, 6, 3, 3).transpose();
    hess.block(9, 0, 3, 3) = hess.block(0, 9, 3, 3).transpose();
    hess.block(6, 3, 3, 3) = hess.block(3, 6, 3, 3).transpose();
    hess.block(9, 3, 3, 3) = hess.block(3, 9, 3, 3).transpose();

    return hess;
}

