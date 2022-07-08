function makePlot(q, imgDir, saveImage, Edges, ElementToEdge, currentTime, h1)

numEdges = numel(Edges) / 2;

% h1 = figure(1);
clf();
hold on
for c=1:numEdges
    %ind = find(ElementToEdge(:) == c);
    %if numel(ind) == 1
    n1 = Edges(c, 1); % node 1 number
    n2 = Edges(c, 2); % node 2 number
    
    x1 = q(3*(n1-1)+1 : 3*n1);
    x2 = q(3*(n2-1)+1 : 3*n2);
    
    plot3( [x1(1), x2(1)], [x1(2), x2(2)], [x1(3), x2(3)], ...
        'b-', 'LineWidth', 2);
    
end

hold off
title( num2str(currentTime, 'Time=%4.2f sec') );
axis equal
box on

xlim([-0.3, 0.3])
ylim([-0.3, 0.3])
zlim([-0.3, 0.0])

view(0, 90)

if saveImage == true
    filename = num2str(currentTime, 't_%4.2f');
    saveas(h1, [imgDir, filename, '.png']);
end

end
