function [Nodes, Edges, ElementToEdge, bendingElements] = generateMesh(l, w, maxMeshSize)
    
minMeshSize = maxMeshSize / 2;

diam = sqrt(w^2 + l^2);
gd = [3; 4; 0; diam/2; 0; -diam/2; 0; diam/2; diam; diam/2];
g = decsg(gd);
model = createpde;
geometryFromEdges(model, g);

FEMesh = generateMesh(model, 'Hmax', maxMeshSize, 'Hmin', minMeshSize, 'GeometricOrder', 'linear');

[~, numNodes] = size(FEMesh.Nodes);
[~, numElements] = size(FEMesh.Elements);
Nodes = FEMesh.Nodes;
% Nodes = [Nodes(1, :); zeros(1, numNodes); -Nodes(2, :)]; 
Nodes = [Nodes; zeros(1, numNodes)];
ElementToNode = FEMesh.Elements;
ElementToEdge = zeros(3, numElements);

numEdges = 0;
Edges = zeros(2, 3*numElements);

for c = 1:numElements
    node1 = min(ElementToNode(:, c));
    node2 = median(ElementToNode(:, c));
    node3 = max(ElementToNode(:, c));

    ind = find(Edges(1,:) == node1 & Edges(2,:) == node2);
    if numel(ind) == 0 % edge not in array yet; add to array
        numEdges = numEdges + 1;
        Edges(1, numEdges) = node1;
        Edges(2, numEdges) = node2;
        edge1 = numEdges;
    elseif numel(ind) == 1 % edge found in array; assign edge1
        edge1 = ind;
    else
        disp('Error in edge 1 location');
    end

    ind = find(Edges(1,:) == node2 & Edges(2,:) == node3);
    if numel(ind) == 0
        numEdges = numEdges + 1;
        Edges(1, numEdges) = node2;
        Edges(2, numEdges) = node3;
        edge2 = numEdges;
    elseif numel(ind) == 1
        edge2 = ind;
    else
        disp('Error in edge 2 location');
    end

    ind = find(Edges(1,:) == node1 & Edges(2,:) == node3);
    if numel(ind) == 0
        numEdges = numEdges + 1;
        Edges(1, numEdges) = node1;
        Edges(2, numEdges) = node3;
        edge3 = numEdges;
    elseif numel(ind) == 1
        edge3 = ind;
    else
        disp('Error in edge 3 location');
    end
    
    ElementToEdge(1, c) = min([edge1 edge2 edge3]);
    ElementToEdge(2, c) = median([edge1 edge2 edge3]);
    ElementToEdge(3, c) = max([edge1 edge2 edge3]);
end

Edges = Edges(:, 1:numEdges); % trim Edges array

numBendElements = 0;
bendingElements = zeros(4, numEdges);

for c = 1:numEdges
    ind = find(ElementToEdge(:) == c);
    if numel(ind) == 1
        % edge is on the boundary
    elseif numel(ind) == 2 
        % edge is not on the boundary (2 elements share this edge)
        
        n0 = Edges(1, c); % one node of this edge
        n1 = Edges(2, c); % other node of this edge

        element1 = ceil(ind(1) / 3); % first element which shares this edge
        element2 = ceil(ind(2) / 3); % second element which shares this edge
        element1_nodes = ElementToNode(:, element1); % get the 3 nodes of first element
        element2_nodes = ElementToNode(:, element2); % get the 3 nodes of second element 
        
        n2 = element1_nodes(element1_nodes~=n0 & element1_nodes~=n1); % node of first element not on the edge
        n3 = element2_nodes(element2_nodes~=n0 & element2_nodes~=n1); % node of second element not on the edge
        
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

Nodes = transpose(Nodes);
Edges = transpose(Edges);
bendingElements = transpose(bendingElements);

end