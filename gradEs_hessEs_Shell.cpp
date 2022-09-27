#include "gradEs_hessEs_Shell.h"

gradEs_hessEs_Shell::gradEs_hessEs_Shell(Vector3d node0, Vector3d node1, double l_k, double ks) {
    Vector3d edge = node1 - node0; 
    double edgeLen = edge.norm();
    Vector3d tangent = edge / edgeLen;
    epsX = edgeLen / l_k - 1;
    dF_unit = ks / l_k * epsX * tangent;
    dF = VectorXd(6);
    dF << - dF_unit(0), - dF_unit(1), - dF_unit(2), dF_unit(0), dF_unit(1), dF_unit(2);
    M = ks / l_k * ((1 / l_k - 1 / edgeLen) * Matrix3d::Identity() + edge * edge.transpose() / pow(edgeLen, 3));
    dJ = MatrixXd(6, 6);
    // dJ << M(0, 0), M(0, 1), M(0, 2), - M(0, 0), - M(0, 1), - M(0, 2),
    //     M(1, 0), M(1, 1), M(1, 2), - M(1, 0), - M(1, 1), - M(1, 2),
    //     M(2, 0), M(2, 1), M(2, 2), - M(2, 0), - M(2, 1), - M(2, 2),
    //     - M(0, 0), - M(0, 1), - M(0, 2), M(0, 0), M(0, 1), M(0, 2),
    //     - M(1, 0), - M(1, 1), - M(1, 2), M(1, 0), M(1, 1), M(1, 2),
    //     - M(2, 0), - M(2, 1), - M(2, 2), M(2, 0), M(2, 1), M(2, 2);
    dJ.block(0, 0, 3, 3) = M;
    dJ.block(0, 3, 3, 3) = - M;
    dJ.block(3, 0, 3, 3) = - M;
    dJ.block(3, 3, 3, 3) = M;
}

gradEs_hessEs_Shell::~gradEs_hessEs_Shell() { 
    ;
}

VectorXd gradEs_hessEs_Shell::getdF() {
    return dF;
}

MatrixXd gradEs_hessEs_Shell::getdJ() {
    return dJ;
}