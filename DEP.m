%% Copyright M. Khalid Jawed (khalidjm@seas.ucla.edu)
% You should use this code at your own risk. Copy and redistribution is not
% permitted. Written permission is required.
%
% Update: 12/13/ 2021

function DEP(h_min, h_max, aspectRatio)
% aspectRatio = Length/Width
% h_min = thickness at the beginning of the simulation
% h_max = thickness at the end of the simulation

%% Parameters
w = 0.20; % width (m)
L = aspectRatio * w; % length
Y = 0.2e6; % Young's modulus

totalTime = 0.05;
numDataPoints = 500; % number of saved configurations

dt = 1e-2; % time step size in seconds

rho = 1200; % density (kg/m^3)

g = [0, 0, -9.81];

visc = 1e-2; % Viscosity

maxIter = 20;

plotStep = 1; % 1000;

saveImage = true;

%% Generate mesh
maxMeshSize = w/15;
[Nodes, Edges, Elements] = generateMesh(L, w, maxMeshSize);

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
    
    iter = 0;
    while err > smallNo && iter < maxIter
        
        % [force, Jac] = getForceJacobian(q);
        getForceJacobian
        
        forceFree = force( freeDOF );
        JacFree = Jac( freeDOF, freeDOF );
        
%         JacFree = sparse(JacFree);
        dq = JacFree \ forceFree; % https://www.pardiso-project.org/
        % JacFree is a square matrix and most of its entries are zero.
        % Pardiso Project will let us take advantage of that.
        
        qFree = qFree - dq;
        
        q( freeDOF ) = qFree;
        
        err = sum( abs( forceFree ) );
        iter = iter + 1;
        % fprintf('Error = %f at iter=%d\n', err, iter);
    end
    
    if iter >= maxIter
        disp('Could not converge');
        return
    end
    
    qDotFree = ( qFree - q0( freeDOF ) ) / dt;
    qDot( freeDOF ) = qDotFree;
    
    maxZ(ct) = max( q(3:3:end) );

    if plotStep > 0 && (mod(ct,plotStep)==0 || ct == numSteps)
        makePlot(q, imgDir, saveImage, Edges, currentTime);
        
        figure(2);
        plot( maxZ (maxZ ~= 0), 'ro' );
        
    end
    
    if mod(ct, numSteps_between_data) == 0
        c_data = c_data + 1;
        xSave(c_data,:) = q;
        hSave(c_data) = h_running;
    end
end

dataName = num2str(aspectRatio, 'data_%4.2f.mat');
save(dataName, 'xSave', 'hSave', 'L', 'w', 'aspectRatio', 'Y', ...
    'totalTime', 'dt', 'rho', 'g', 'visc');

end
