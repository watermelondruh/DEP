#ifndef GRADES_HESSES_SHELL_H
#define GRADES_HESSES_SHELL_H

#include <Eigen/Dense>
using namespace Eigen;

class gradEs_hessEs_Shell
{
public:
	gradEs_hessEs_Shell(Vector3d node0, Vector3d node1, double l_k, double ks);
	~gradEs_hessEs_Shell();
    VectorXd getdF();
    MatrixXd getdJ();

private:
    double l_k;
    double ks;
    double epsX;
    
    Vector3d dF_unit;
    Vector3d node0;
    Vector3d node1;
    VectorXd dF;
    MatrixXd dJ;
    Matrix3d M;
};

#endif
