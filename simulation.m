clear all;
close all;
clc;

l = 0.2; % length (m)
w = 0.2; % width (m)
h_min = 0.0015; % min thickness (m)
h_max = 0.0015; % max thickness (m)
Y = 0.18e6; % Young's modulus (Pa)
rho = 1020; % density (kg/m^3)
visc = 1e-2; % viscosity (kg/(ms))
g = [0, 0, -9.81]; % gravity (m/s^2)

totalTime = 10; % (s)
dt = 0.01; % (s)
numDataPoints = 500;
maxIter = 20;
plotStep = 1; 

% numElements = 2;
% numNodes = 4;
% numEdges = 5;
% Elements = [2, 4, 1, 3];
% Edges = [1,2; 2,4; 1,4; 2,3; 3,4];
% Nodes = [0,0,0; 0.1,0,0; 0.1,0.1,0; 0,0.1,0]; % [0,0,0; 0.1,0,-0.1; 0,0,-0.2; -0.1,0,-0.1];

maxMeshSize = w / 8;
[Nodes, Edges, ElementsToEdge, Elements] = generateMesh(l, w, maxMeshSize);
[numNodes, ~] = size(Nodes);
[numEdges, ~] = size(Edges);
[numElements, ~] = size(Elements);

% perturb slightly
R = rotx(-80);
Nodes = Nodes * R';

% rotate by 90 degrees
% Nodes = [Nodes(:, 1), Nodes(:, 3), -Nodes(:, 2)];

% perturb slightly
% diam = sqrt(w^2 + l^2);
% for c = 1 : numNodes
%     Nodes(c, 2) = Nodes(c, 2) + 0.05 * abs(abs(Nodes(c, 3)) - diam / 2);
%     % Nodes(c, 2) = Nodes(c, 2) + 0.05 * abs(Nodes(c, 3));
%     % Nodes(c, 2) = 0.15 * sin(Nodes(c, 1) * 2 * pi / diam);
% end

q = zeros(3*numNodes, 1);
for c = 1:numNodes
    q(3*c-2:3*c) = Nodes(c, :);
end

MassVector = zeros(3*numNodes, 1);
totalM = l * w * h_max * rho;
dm = totalM / numNodes;
MassVector(:) = dm;
MassMatrix = diag(MassVector);

% weight vector
W = zeros(3*numNodes, 1);
for c = 1:numNodes
    W(3*c-2:3*c) = dm * g;
end

% reference length of each edge
refLen = zeros(numEdges, 1);
for c = 1:numEdges
    n0 = Edges(c, 1);
    n1 = Edges(c, 2);
    refLen(c) = norm(Nodes(n0,:) - Nodes(n1,:));
end

viscLen = mean(refLen);

ks = 1/2 * sqrt(3) * Y * h_max * refLen.^2; % stretching stiffness: varies (w.r.t. refLen)
kb = 2/sqrt(3) * Y *(h_max^3)/12; % bending stiffness: uniform
tolerance = kb / l * 1e-2; % tolerance (divide by L instead of L^2)

freeDOF = [7:9 13:3*numNodes];
% freeDOF = [4:3*numNodes];

qDot = q * 0; % velocity
qDot0 = qDot; % initial velocity

numSteps = floor(totalTime / dt);
numSteps_between_data = round(numSteps/numDataPoints);

maxZ = zeros(numSteps, 1);

xSave = zeros(numDataPoints+1, 3*numNodes);
hSave = zeros(numDataPoints+1, 1);
c_data = 1;
xSave(c_data,:) = q;
hSave(c_data) = h_max;

k = 1;

holdTime = 50;

for ct = 1:numSteps
    
    if ct==holdTime/5
        freeDOF = [(7:3*numNodes)]';
    elseif ct>holdTime
        freeDOF = (4:3*numNodes)';
    end

    currentTime = dt * ct;
    
    q0 = q; % initial position
    qDot0 = qDot; % initial velocity
    qFree = q0(freeDOF);

    err = 10 * tolerance;
    iter = 0;

    while err > tolerance && iter < maxIter
        getForceJacobian
        
        forceFree = force(freeDOF);
        JacFree = Jac(freeDOF, freeDOF);
        
        % Newton's Method:
        % - scalar form: x^{k+1} <- x^k - f(x^k) / f'(x^k)
        % - in matrix form, f(x^k) \ f'(x^k) ==> Jacobian \ Force

        % JacFree = sparse(JacFree);
        dq = JacFree \ forceFree; % https://www.pardiso-project.org/
        % JacFree is a square matrix and most of its entries are zero.
        % Pardiso Project will let us take advantage of that.
        
        qFree = qFree - dq;
        q(freeDOF) = qFree;
        
        err = sum(abs(forceFree));
        iter = iter + 1;

    end
    
    if iter >= maxIter
        disp('Could not converge');
        return
    end
    
    qDotFree = (qFree - q0(freeDOF)) / dt;
    qDot(freeDOF) = qDotFree;
    
    maxZ(ct) = max(q(3:3:end));
    
    filename = sprintf("l=%.01f_w=%.01f_h=%.04f_hold=%i_dt=%.02f.gif", l, w, h_max, holdTime, dt);
    imgdir = './';
    del = 0.01;

    if plotStep > 0 && (mod(ct,plotStep)==0 || ct == numSteps)
        h1 = figure(1);
        makePlot(q, imgdir, false, Edges, ElementsToEdge, currentTime, h1);
        xlabel('x');
        ylabel('y');
        zlabel('z');

        frame = getframe(h1); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256); 
        
        % write .gif 
        if k == 1
          imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',del); 
        else 
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',del); 
        end  
        
    end
    
    k = k + 1;
    
    if mod(ct, numSteps_between_data) == 0
        c_data = c_data + 1;
        xSave(c_data,:) = q;
        hSave(c_data) = h_running;
    end

end
