#include <bits/stdc++.h>

int main() {
    double w = 0.20; // width (m)
    double L = aspectRatio * w; // length
    double Y = 0.2e6; // Young's modulus

    double totalTime = 0.05;
    int numDataPoints = 500; // number of saved configurations

    double dt = 1e-2; // time step size in seconds

    double rho = 1200; // density (kg/m^3)

    int[] g = [0, 0, -9.81];
    int visc = 1e-2;
    int maxIter = 20;
    int plotStep = 1; % 1000
    bool saveImage = true;

    double maxMeshSize = w/15;


%%
% Nodes = load('Nodes.txt', 'w');
[numNodes,~] = size(Nodes);

% Edges = load('Edges.txt', 'w');
[numEdges,~] = size(Edges);

% Elements = load('Elements.txt', 'w');
[numElements,~] = size(Elements);

% DOF vector
q = zeros(3*numNodes, 1);
for c=1:numNodes
    q( 3*(c-1)+1: 3*c ) = Nodes(c,:);
end

% Mass vector
Mass_vec = zeros( 3*numNodes, 1);
totalM = L * w * h_max * rho;
dm = totalM / numNodes;
Mass_vec(:) = dm;
MassMatrix = diag(Mass_vec);

% Weight vector
W = zeros( 3*numNodes, 1);
for c=1:numNodes
    W( 3*(c-1)+1: 3*c ) = dm * g;
end

% Reference length of each edge
refLen = zeros(numEdges, 1);
for c=1:numEdges
    n0 = Edges(c,1);
    n1 = Edges(c,2);
    refLen(c) = norm( Nodes(n0,:) - Nodes(n1,:) );
end

viscLen = mean(refLen);

%% Bending and stretching stiffness
ks = 1/2 * sqrt(3) * Y * h_max * refLen.^2;
kb = 2/sqrt(3) * Y *(h_max^3)/12;

%% Tolerance
smallNo = kb / L * 1e-2;

%% Fixed and Free DOFs
fixedDOF = 1:3;
dummy = ones(3*numNodes, 1);
dummy(fixedDOF) = 0;
freeDOF = find( dummy == 1 );

%%
imgDir = [date, '/'];
if saveImage == true
    mkdir(imgDir);
end

%% Time marching
qDot = q * 0; % velocity
qDot0 = qDot; % Initial velocity

time_arr = 0:dt:totalTime;
numSteps = numel(time_arr);
numSteps_between_data = round(numSteps/numDataPoints);

maxZ = zeros(numSteps,1);

xSave = zeros(numDataPoints+1, 3*numNodes);
hSave = zeros(numDataPoints+1, 1);
c_data = 1;
xSave(c_data,:) = q;
hSave(c_data) = h_max;

for ct = 2:numSteps
    
    currentTime = time_arr(ct);
    
    if mod(ct,100) == 0
        fprintf('Time=%f aspectRatio=%f\n', currentTime, aspectRatio);
    end
    
    q0 = q; % Initial position
    qDot0 = qDot; % Initial velocity
    
    err = 10 * smallNo;
    
    qFree = q0 ( freeDOF );
}