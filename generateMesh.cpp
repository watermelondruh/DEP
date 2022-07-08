#include <bits/stdc++.h>

void generateMesh(double L, double w, doube maxMeshSize) {
    // return nodes, edges, elements
    double l = 0.1;
    double w = 0.03;
    double maxMeshSize = w/1;
    double minMeshSize = maxMeshSize/2;

    double[] gd = [3; 4; 0; l; l; 0; 0; 0; w; w];
    double g = decsg(gd);
    Model model = createpde; // todo: make model
    geometryFromEdges(model,g);

FEMesh = generateMesh(model,'Hmax', maxMeshSize, 'Hmin', ...
    minMeshSize, 'GeometricOrder', 'linear');

%%
ElementToNode = FEMesh.Elements;
Nodes = FEMesh.Nodes;
[~, numNodes] = size(FEMesh.Nodes);
[~, numElements] = size(FEMesh.Elements);

% Convert 2xN matrix to 3xN matrix
Nodes = [Nodes; zeros(1, numNodes)];

ElementToEdge = zeros(3, numElements);
Edges = zeros( 2, 3 * numElements); % large enough container

numEdges = 0;
for cEl = 1 : numElements
    
    node1 = min( ElementToNode(:, cEl) );
    node2 = median( ElementToNode(:, cEl) );
    node3 = max( ElementToNode(:, cEl) );
    
    % Go over the three edges and add to "Edges" matrix if necessary
    ind = find( Edges(1,:) == node1 & Edges(2,:) == node2 );
    if numel(ind) == 0
        numEdges = numEdges + 1;
        Edges( 1, numEdges ) = node1;
        Edges( 2, numEdges ) = node2;
        
        edge1_no = numEdges;
    elseif numel(ind) == 1
        edge1_no = ind;
    else
        disp('Error in edge 1 location');
    end
    
    ind = find( Edges(1,:) == node2 & Edges(2,:) == node3 );
    if numel(ind) == 0
        numEdges = numEdges + 1;
        Edges( 1, numEdges ) = node2;
        Edges( 2, numEdges ) = node3;
        
        edge2_no = numEdges;
    elseif numel(ind) == 1
        edge2_no = ind;
    else
        disp('Error in edge 2 location');
    end
    
    ind = find( Edges(1,:) == node1 & Edges(2,:) == node3 );
    if numel(ind) == 0
        numEdges = numEdges + 1;
        Edges( 1, numEdges ) = node1;
        Edges( 2, numEdges ) = node3;
        
        edge3_no = numEdges;
    elseif numel(ind) == 1
        edge3_no = ind;
    else
        disp('Error in edge 3 location');
    end
    // element: finite element analysis
    // node: corners of an element
    ElementToEdge( 1, cEl ) = min( [edge1_no edge2_no edge3_no] );
    ElementToEdge( 2, cEl ) = median( [edge1_no edge2_no edge3_no] );
    ElementToEdge( 3, cEl ) = max( [edge1_no edge2_no edge3_no] );
end

Edges = Edges( :, 1:numEdges);

%% Construct Bending Elements
bendingElements = zeros(4, numEdges);

% Adapted from https://github.com/shift09/plates-shells/blob/master/src/bending.cpp
% //         x2
% //         /\
% //        /  \
% //     e1/    \e3
% //      /  t0  \
% //     /        \
% //    /    e0    \
% //  x0------------x1
% //    \          /
% //     \   t1   /
% //      \      /
% //     e2\    /e4
% //        \  /
% //         \/
% //         x3

numBendElements = 0;
for c=1:numEdges
    
    % Find the two triangular elements that share this edge
    ind = find( ElementToEdge(:) == c );
    if numel(ind) == 1
        % fprintf('Edge %d is on the boundary\n', c);
    elseif numel(ind) == 2
        n0 = Edges(1, c);
        n1 = Edges(2, c);
        
        El_1 = ceil( ind(1) / 3 ); % Triangle 1 that includes edge no. c
        El_2 = ceil( ind(2) / 3 ); % Triangle 2 that includes edge no. c
        
        El_1_nodes = ElementToNode(:, El_1);
        El_2_nodes = ElementToNode(:, El_2);
        
        n2 = El_1_nodes( El_1_nodes~=n0 & El_1_nodes~=n1 );
        n3 = El_2_nodes( El_2_nodes~=n0 & El_2_nodes~=n1 );
        
        % Store in the container
        numBendElements = numBendElements + 1;
        bendingElements(1, numBendElements) = n0;
        bendingElements(2, numBendElements) = n1;
        bendingElements(3, numBendElements) = n2;
        bendingElements(4, numBendElements) = n3;
        
    else
        disp('Error in finding bending elements');
    end
end

bendingElements = bendingElements(:, 1:numBendElements);

figure(2);
clf();
pdemesh(FEMesh.Nodes, FEMesh.Elements,'EdgeColor','k');
hold on
for c=1:numEdges
    n1 = Edges(1, c); % node 1 number
    n2 = Edges(2, c); % node 2 number
    
    plot( [Nodes(1,n1), Nodes(1,n2)], [Nodes(2,n1), Nodes(2,n2)], ...
        'b--', 'LineWidth', 3);
    
end

for c=1:numBendElements
    x0 = Nodes(:, bendingElements(1,c));
    x1 = Nodes(:, bendingElements(2,c));
    x2 = Nodes(:, bendingElements(3,c));
    x3 = Nodes(:, bendingElements(4,c));
    
    plot( [x0(1), x1(1)], [x0(2), x1(2)], 'r^-', 'LineWidth', 2);
    plot( [x0(1), x2(1)], [x0(2), x2(2)], 'r--');
    plot( [x0(1), x3(1)], [x0(2), x3(2)], 'r--');
    plot( [x1(1), x2(1)], [x1(2), x2(2)], 'r--');
    plot( [x1(1), x3(1)], [x1(2), x3(2)], 'r--');
end

hold off

%% Save
saveData = 0;
if saveData == 1
    fid = fopen('Nodes.txt', 'w');
    for c=1:numNodes
        for d=1:3
            fprintf(fid, '%d ', Nodes(d,c) );
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
    
    fid = fopen('Edges.txt', 'w');
    for c=1:numEdges
        for d=1:2
            fprintf(fid, '%d ', Edges(d,c) );
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
    
    fid = fopen('Elements.txt', 'w');
    for c=1:numBendElements
        for d=1:4
            fprintf(fid, '%d ', bendingElements(d,c) );
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
end

Nodes = transpose(Nodes);
Edges = transpose(Edges);
bendingElements = transpose(bendingElements);

end