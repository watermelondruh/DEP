#ifndef GRADEB_HESSEB_SHELL_H
#define GRADEB_HESSEB_SHELL_H

#include "theta.h"
#include <Eigen/Dense>
using namespace Eigen;

class gradEb_hessEb_Shell
{
public:
	gradEb_hessEb_Shell(Vector3d x0, Vector3d x1, Vector3d x2, Vector3d x3, double thetabar, double kb);
	~gradEb_hessEb_Shell();
    VectorXd getdF();
    MatrixXd getdJ();

private:
    VectorXd dF;
    MatrixXd dJ;
    double theta;
    VectorXd grad;
    MatrixXd hess;
};

#endif
