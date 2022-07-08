/* NOTE: last modified in February */

#include <bits/stdc++.h>
#include "utils.h"

#define sqrt3 1.73205080757

using namespace std;

float h_min = 0.2;
float h_max = 0.2;
float aspectRatio = 0.2;
float w = 0.2; // width (m)
float L = aspectRatio * w; // length (m)
float Y = 0.2e6; // Young's modulus
float dt = 1e-2; // time step size in seconds
float rho = 1200; // density (kg/m^3)
float g[3] = {0, 0, -9.81};
float visc = 1e-2; // viscosity

float totalTime = 10;
int numDataPoints = 500; // number of saved configurations
int maxIter = 20;
int plotStep = 1; // 1000
bool saveImage = false;

string filename = "testAnimated.gif";
float del = 0.01;

int main() {

    /* --- GENERATE MESH --- */
    float maxMeshSize = w / 15;
    // [Nodes, Edges, Elements] = generateMesh(L, w, maxMeshSize);
    // TODO: modify vectors within functions
    vector<vector<vector<float>>> NodesEdgesElements = generateMesh(L, w, maxMeshSize);
    // Nodes = load('Nodes.txt', 'w');
    vector<float> nodes = NodesEdgesElements[0];
    int numNodes = Nodes[0].size();
    // Edges = load('Edges.txt', 'w');
    vector<float> edges = NodesEdgesElements[1];
    int numEdges = Edges[0].size();
    // Elements = load('Elements.txt', 'w');
    vector<float> elements = NodesEdgesElements[2];
    int numElements = Elements[0].size();


    /* --- DOF VECTOR --- */
    vector<float> q(3 * numNodes, 0.0);
    vector<float> q0(3 * numNodes, 0.0);
    // vector<vector<float>> q(3 * numNodes, vector<float>(1, 0.0));
    for (int c = 0; c < numNodes; c++) {
        for (int i = 0; i < 3; i++) {
            q[3 * c + i] = Nodes[c][i];
            // q[3 * c + i][0] = Nodes[c][i];
        }
    }
    vector<float> qFree(3 * numNodes - 3, 0.0);


    /* --- MASS VECTOR --- */
    vector<float> MassVec(3 * numNodes, 0.0);
    // vector<vector<float>> MassVec(3 * numNodes, vector<float>(1, 0.0));
    vector<vector<float>> MassMatrix(3 * numNodes, vector<float>(3 * numNodes, 0.0));
    float totalM = L * w * h_max * rho;
    float dm = totalM / numNodes;
    for (int c = 0; c < 3 * numNodes; c++) {
        MassVec[c] = dm;
        // MassVec[c][0] = dm;
        MassMatrix[c][c] = dm;
    }


    /* --- WEIGHT VECTOR --- */
    vector<float> W(3 * numNodes, 0.0);
    // vector<vector<float>> W(3 * numNodes, vector<float(1, 0.0));
    for (int c = 0; c < numNodes; c++) {
        for (int i = 0; i < 3; i++) {
            W[3 * c + i] = dm * g[i];
            // W[3 * c + i][0] = dm * g[i];
        }
    }


    /* --- REFERENCE LENGTH OF EACH EDGE --- */
    vector<float> refLen(numEdges, 0.0);
    // vector<vector<float>> refLen(numEdges, vector<float>(1, 0.0));
    for (int c = 0; c < numEdges; c++) {
        int n0 = (int) Edges[c, 0];
        int n1 = (int) Edges[c, 1];
        refLen[c] = norm(Nodes[n0] - Nodes[n1]);
        // refLen[c][0] = norm(Nodes[n0] - Nodes[n1]);
    }

    float viscLen = mean(refLen);

    /* --- BENDING AND STRETCHING STIFFNESS --- */
    vector<float> ks(numEdges, 0.0);
    // vector<vector<float>> ks(numEdges, vector<float>(1, 0.0));
    float temp1 = sqrt3 * Y * h_max / 2;
    for (int c = 0; c < numEdges; c++) {
        ks[c] = temp1 * pow(refLen[c], 2);
        // ks[c][0] = tmep1 * pow(refLen[c][0], 2)
    }
    float kb = 2 * sqrt3 * Y * pow(h_max, 3) / 12;


    /* --- TOLERANCE --- */
    float smallNo = kb / L * 1e-2;


    /* -- FIXED AND FREE DOFS -- */
    vector<int> dummy(3 * numNodes, 1);
    // vector<vector<int>> dummy(3 * numNodes, vector<int>(1, 1));
    for (int i = 0; i < 3; i++) {
        dummy[i] = 0;
        //dummy[i][0] = 0;
    }
    vector<int> freeDOF(3 * numNodes - 3, 0);
    for (int i = 0; i < 3 * numNodes - 3; i++) {
        freeDOF[i] = i + 3;
    }

    
    /* --- IMAGE DIRECTORY --- */
    /* imgDir = [date, '/'];
        if saveImage == true
            mkdir(imgDir);
        end */


    /* --- TIME MARCHING --- */
    vector<float> qDot(3 * numNodes, 0.0); // velocity
    vector<float> qDot0(3 * numNodes, 0.0); // initial velocity
    // vector<vector<float>> qDot(3 * numNodes, vector<float>(1, 0.0)); 
    // vector<vector<float>> qDot0(3 * numNodes, vector<float>(1, 0.0)); 

    int numSteps = (int) (totalTime / dt);
    int numSteps_between_data = (int) roundf(1.0 * numSteps / numDataPoints);

    vector<float> maxZ(numSteps, 0.0);
    // vector<vector<float>> maxZ(numSteps, vector<float>(1, 0.0));

    int c_data = 0;
    vector<vector<float>> xSave(numDataPoints + 1, vector<float>(3 * numNodes, 0.0));
    vector<float> hSave(numDataPoints + 1, 0.0);
    // vector<vector<float>> hSave(numDataPoints + 1, vector<float>(1, 0.0));
    
    for (int c = c_data; c < 3 * numNodes; c++) {
        xSave[c][0] = q[c][0];
        hSave[c] =  h_max;
        // hSave[c][0] = h_max;
    }

    
    vector<float> force;
    vector<vector<float>> Jac;
    vector<float> forceFree;
    vector<vector<float>> JacFree;

    int k = 1;

    for (int ct = 1; ct < numSteps; ct++) {
        float currentTime = ct * dt;
        if (ct % 100 == 0) {
            cout << 'Time = ' << currentTime << ', aspectRatio = ' << aspectRatio << endl;
        }
        
        for (int i = 0; i < 3 * numNodes; i++) {
            q0[i] = q[i]; // initial position
            qDot0[i] = qDot[i]; // initial velocity
        }

        float err = 10 * smallNo;
        
        for (int i = 3; i < 3 * numNodes; i++) {
            qFree[i - 3] = q0[i];
        }

        int iter = 0;
        while (err > smallNo && iter < maxIter) {
            // TODO: modify vector directly in getForceJacobian
            getForceJacobian(q, &force, &Jac);
            
            // TODO: correct dimensions
            for (int i = 3; i < force.size(); i++) {
                forceFree[i] = force[i];
                for (int j = 3; j < force[0].size(); j++) {
                    JacFree[i][j] = Jac[i][j];
                }
            }

            // JacFree = sparse(JacFree);
            // dq = JacFree \ forceFree; // https://www.pardiso-project.org/
            /* "JacFree is a square matrix and most of its entries are zero.
                Pardiso Project will let us take advantage of that." */

            for (int i = 3; i < 3 * numNodes; i++) {
                qFree[i - 3] -= dq;
                q[i] = qFree[i - 3];
            }

            for (int i = 0; i < forceFree.size(); i++) {
                err += abs(forceFree[i]);
            }

            iter += 1;
            
            cout << "iter = " << iter << ", error = " << err << endl;
        }

        if (iter >= maxIter) {
            cout << "Could not converge." << endl;
            return;
        }

        for (int i = 3; i < 3 * numNodes; i++) {
            qDotFree[i - 3] = (qFree[i - 3] - q0[i]) / dt;
            qDot[i] = qDotFree[i - 3];
        }

        // TODO: stride = 3
        // maxZ(ct) = max( q(3:3:end) );
        maxZ[ct] = *max_element(q.begin() + 3; q.end());

        if ((plotStep > 0 && ct % plotStep == 0) || ct == numSteps) {
            // TODO: plotting in C++
            /* 
            h1 = figure(1);
            makePlot(q, imgDir, saveImage, Edges, currentTime, h1);

            frame = getframe(h1); 
            im = frame2im(frame); 
            [imind,cm] = rgb2ind(im,256); 
            
            % Write to the GIF File 
            if k == 1 
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',del); 
            else 
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',del); 
            end  
            
    %         figure(2);
    %         plot( maxZ (maxZ ~= 0), 'ro' );
            */
        }

        k += 1;

        if (ct % numSteps_between_data == 0) {
            c_data += 1;
            for (int c = c_data; c < numDataPoints + 1; c++) {
                for (int i = 0; i < 3 * numNodes; i++) {
                    xSave[c][i] = q[i];
                }
                hSave[c] = h_running;
            }
        }
    }

    // TODO: save data
    /* 
    dataName = num2str(aspectRatio, 'data_%4.2f.mat');
    save(dataName, 'xSave', 'hSave', 'L', 'w', 'aspectRatio', 'Y', ...
    'totalTime', 'dt', 'rho', 'g', 'visc'); 
    */
}

