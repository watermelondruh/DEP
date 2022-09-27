#ifndef THETA_H
#define THETA_H

#include <Eigen/Dense>
using namespace Eigen;

class Theta {
public:
	Theta(Vector3d x0, Vector3d x1, Vector3d x2, Vector3d x3);
	~Theta();

    double getTheta();
    VectorXd gradTheta();
    MatrixXd hessTheta();

private:
    Vector3d m_e0, m_e1, m_e2, m_e3, m_e4;
    Vector3d m_nn1, m_nn2;

    double m_cosA1, m_cosA2, m_cosA3, m_cosA4;
    double m_sinA1, m_sinA2, m_sinA3, m_sinA4;
    double m_h1, m_h2, m_h3, m_h4, m_h01, m_h02;
};

#endif
