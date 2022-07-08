function grad_hess_elastic(q)
    N = numel(q) / 3;
    F = zeros(3*N, 1);
    G = zeros(3*N, 3*N);

    ind = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];

    numEdges = N / 3;
    for k = 1 : numEdges
        x0 = Edges(1, k)
        x1 = Edges(2, k)
        ind = ?
        [dF, dJ] = gradEs_hessEs_shell(x0, x1);
        F(ind) = f(ind) + dF;
        J(ind) = J(ind, ind) + dJ;
    end

    numHinges = ?
    for k = 1 : numHinges
        x0 = 
    