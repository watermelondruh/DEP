#include "gradEb_hessEb_Shell.h"

gradEb_hessEb_Shell::gradEb_hessEb_Shell(Vector3d x0, Vector3d x1, Vector3d x2, Vector3d x3, double thetaBar, double kb) {
    dF = VectorXd::Zero(12);
    dJ = MatrixXd::Zero(12, 12);
    Theta Theta(x0, x1, x2, x3);
    theta = Theta.getTheta();
    grad = Theta.gradTheta();
    dF = 0.5 * kb * (2 * (theta - thetaBar) * grad);
    hess = Theta.hessTheta();
    dJ = 0.5 * kb * (2 * grad * grad.transpose() + 2 * (theta - thetaBar) * hess);
}

gradEb_hessEb_Shell::~gradEb_hessEb_Shell() { 
    ;
}

VectorXd gradEb_hessEb_Shell::getdF() {
    return dF;
}

MatrixXd gradEb_hessEb_Shell::getdJ() {
    return dJ;
}