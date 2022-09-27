#ifndef SIMULATION_H
#define SIMULATION_H

#include "gradEb_hessEb_Shell.h"
#include "gradEs_hessEs_Shell.h"
#include "mfem/mfem.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <iostream>
#include <fstream>

using namespace std;
using namespace mfem;
using namespace Eigen;

class Simulation
{
public:
    Simulation(double l, double w, double h,
               double Y, double rho, double visc, 
               double totalTime, double dt);
    ~Simulation();
    void generateMesh();
	void getForceJacobian();
    void run();

private:
    Matrix3Xd Nodes;
    Matrix2Xi Edges;
    Matrix3Xi ElementToNode;
    Matrix3Xi ElementToEdge;
    Matrix4Xi BendingElements;
    
    VectorXd q, q0, qDot, qDot0, dq;
    VectorXd refLen;
    VectorXd MassVector;
    MatrixXd MassMatrix;
    VectorXd Fs, Fb, Fv, Force, ForceFree;
    MatrixXd Js, Jb, Jv, Jac, JacFree;
    VectorXd ks, ks_running;
    VectorXd W, W_running;

    int numNodes, numElements, numEdges;
    double kb, kb_running;
    double visc, viscLen;
    double currentTime, totalTime, dt;
    double h_running, h_min, h_max;
    double l, w;
    double Y, rho;
};

#endif
