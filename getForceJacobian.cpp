#include <bits/stdc++.h>
#include "utils.h"

int h_running() {
    if (currentTime <= totalTime / 2)
        return h_min + (h_max - h_min) * currentTime / (totalTime / 2);
    else
        return h_max + (h_min - h_max) * (currentTime - totalTime / 2 ) / (totalTime / 2);
}

int getForceJacobian(int q) {
    /*
% global kb ks
% global q0 qDot0
% global dt
% global Edges Elements
% global refLen
% global MassMatrix Mass_vec
% global W
% global visc viscLen
    */

    int ndof = q.size();
    int numEdges = Edges[0].size();

    // Compute stretching force
    vector<int> Fs(ndof, 0);
    vector<vector<int>> Js(ndof, vector<int>(ndof, 0));

    for (int c = 0; c < numEdges; c++) {
        int n0 = Edges[c][0], n1 = Edges[c][1];
        // node0 = q(3*(n0-1)+1:3*n0);
        // node1 = q(3*(n1-1)+1:3*n1);
        vector<int> node0 = vector<int>(q.begin() + 3 * (n0 - 1), q.begin() + 3 * n0);
        vector<int> node1 = vector<int>(q.begin() + 3 * (n1 - 1), q.begin() + 3 * n1) ;
    
        // ind = [3*(n0-1)+1:3*n0 3*(n1-1)+1:3*n1];
        node0 = transpose(node0);
        node1 = transpose(node1);

        /*
        [dFs, dJs] = gradEs_hessEs_Shell(node0, node1, refLen(c), ks_running(c));
    
        Fs( ind ) = Fs ( ind ) + dFs;
        Js( ind, ind ) = Js ( ind, ind ) + dJs;  
        */
    }
    
    
      
end

% Compute bending force
[numElements,~] = size(Elements);
Fb = zeros(ndof, 1);
Jb = zeros(ndof, ndof);
for c=1:numElements
    n0 = Elements(c,1);
    n1 = Elements(c,2);
    n2 = Elements(c,3);
    n3 = Elements(c,4);
    
    node0 = q(3*(n0-1)+1:3*n0);
    node1 = q(3*(n1-1)+1:3*n1);
    node2 = q(3*(n2-1)+1:3*n2);
    node3 = q(3*(n3-1)+1:3*n3);
    
    ind = [3*(n0-1)+1:3*n0 3*(n1-1)+1:3*n1 3*(n2-1)+1:3*n2 3*(n3-1)+1:3*n3];
    
    thetaBar = 0;
    [dFb, dJb] = gradEb_hessEb_Shell(node0, node1, node2, node3, thetaBar, kb_running);
    Fb( ind ) = Fb ( ind ) + dFb;
    Jb( ind, ind ) = Jb ( ind, ind ) + dJb;    
end

% Compute viscous force
%viscous drag
Fv = -visc * (q-q0)/dt * viscLen;
Jv = -visc * viscLen / dt * eye(ndof,ndof);

force = Mass_vec/dt .* ( (q-q0)/dt - qDot0 ) + Fs + Fb - Fv - W_running;

Jac = MassMatrix/dt^2 + Js + Jb - Jv;

% save('data2.mat');
}