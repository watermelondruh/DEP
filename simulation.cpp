#include "simulation.h"

Simulation::Simulation(double l, double w, double h,
                       double Y, double rho, double visc,
                       double totalTime, double dt)
{
    l = l;
    w = w;
    h_min = h;
    h_max = h;
    Y = Y;
    rho = rho;
    visc = visc;
    totalTime = totalTime;
    dt = dt;
}

Simulation::~Simulation()
{
    ;
}

void Simulation::generateMesh()
{
    const char *mesh_file = "meshes/square.vtk";
    Mesh *mesh = new Mesh(mesh_file, 1, 1);

    numNodes = mesh->GetNV();
    Nodes = Matrix3Xd(3, numNodes);
    for (int i = 0; i < numNodes; i++)
    {
        double *coordinates = mesh->GetVertex(i);
        Nodes(0, i) = coordinates[0];
        Nodes(1, i) = coordinates[1];
        Nodes(2, i) = coordinates[2];
    }

    numElements = mesh->GetNE();
    ElementToNode = Matrix3Xi(3, numElements);
    for (int i = 0; i < numElements; i++)
    {
        Element *element = mesh->GetElement(i);
        mfem::Array<int> vertices;
        element->SetAttribute(1);
        element->GetVertices(vertices);
        ElementToNode(0, i) = vertices[0];
        ElementToNode(1, i) = vertices[1];
        ElementToNode(2, i) = vertices[2];
    }

    numEdges = 0;
    Edges = Matrix2Xi::Zero(2, 3 * numElements);
    ElementToEdge = Matrix3Xi(3, numElements);
    for (int i = 0; i < numElements; i++)
    {
        int node1 = (int) min(min(ElementToNode(0, i), ElementToNode(1, i)), ElementToNode(2, i));
        int node2 = (int) max(min(ElementToNode(0, i), ElementToNode(1, i)), min(max(ElementToNode(0, i), ElementToNode(1, i)), ElementToNode(2, i)));
        int node3 = (int) max(max(ElementToNode(0, i), ElementToNode(1, i)), ElementToNode(2, i));
        int count = 0;
        int edge1 = -1, edge2 = -1, edge3 = -1;
        int j = 0;
        for (j = 0; j < 3 * numElements; j++)
        {
            if (Edges(0, j) == node1 && Edges(1, j) == node2)
            {
                count++;
                edge1 = j;
            }
        }
        if (count == 0) // did not find current edge in array
        { 
            Edges(0, numEdges) = node1;
            Edges(1, numEdges) = node2;
            edge1 = numEdges;
            numEdges++;
        }
        else if (count == 1) // found current edge in array
        { 
            count = 0;
        } 
        else // found multiple instances of current edge in array
        { 
            cout << "Error in edge 1 location." << endl;
            return;
        }
        for (j = 0; j < 3 * numElements; j++)
        {
            if (Edges(0, j) == node2 && Edges(1, j) == node3)
            {
                count++;
                edge2 = j;
            }
        }
        if (count == 0) // did not find current edge in array
        { 
            Edges(0, numEdges) = node2;
            Edges(1, numEdges) = node3;
            edge2 = numEdges;
            numEdges++;
        }
        else if (count == 1) // found current edge in array
        { 
            count = 0;
        }
        else // found multiple instances of current edge in array
        { 
            cout << "Error in edge 2 location." << endl;
            return;
        }
        for (j = 0; j < 3 * numElements; j++)
        {
            if (Edges(0, j) == node1 && Edges(1, j) == node3)
            {
                count++;
                edge3 = j;
            }
        }
        if (count == 0) // did not find current edge in array
        { 
            Edges(0, numEdges) = node1;
            Edges(1, numEdges) = node3;
            edge3 = numEdges;
            numEdges++;
        }
        else if (count == 1) // found current edge in array
        { 
            count = 0;
        }
        else // found multiple instances of current edge in array
        { 
            cout << "Error in edge 3 location." << endl;
            return;
        }
        ElementToEdge(0, i) = (int) min(min(edge1, edge2), edge3);
        ElementToEdge(1, i) = (int) max(min(edge1, edge2), min(max(edge1, edge2), edge3));
        ElementToEdge(2, i) = (int) max(max(edge1, edge2), edge3);
    }
    // Matrix2Xi EdgesN(2, numEdges);
    // for (int i = 0; i < numEdges; i++)
    // {
    //     EdgesN(0, i) = Edges(0, i);
    //     EdgesN(1, i) = Edges(1, i);
    //     cout << EdgesN(0, i) << " " << EdgesN(1, i) << endl;
    // }
    cout << "Number of edges counted: " << numEdges << endl;
    cout << "Number of edges according to MFEM: " << mesh->GetNEdges() << endl;

    int numBendingElements = 0;
    BendingElements = Matrix4Xi(4, numEdges);
    for (int i = 0; i < numEdges; i++)
    {
        vector<int> indices;
        for (int j = 0; j < numElements; j++)
        {
            if (ElementToEdge(0, j) == i || ElementToEdge(1, j) == i || ElementToEdge(2, j) == i)
            {
                indices.push_back(j);
            }
        }
        if (indices.size() == 1) // edge is on the boundary
        { 
            ;
        }
        else if (indices.size() == 2) // edge is not on the boundary; two elements share this edge
        { 
            int node0 = Edges(0, i);
            int node1 = Edges(1, i);
            int element1 = (int) ceil(indices[0] / 3);
            int element2 = (int) ceil(indices[1] / 3);
            Vector3i element1_nodes = ElementToNode.col(element1);
            Vector3i element2_nodes = ElementToNode.col(element2);
            int node2 = -1, node3 = -1;
            for (int j = 0; j < 3; j++)
            {
                if (element1_nodes(j) != node0 && element1_nodes(j) != node1)
                {
                    node2 = element1_nodes(j);
                }
                if (element2_nodes(j) != node0 && element2_nodes(j) != node1)
                {
                    node3 = element2_nodes(j);
                }
            }
            BendingElements(0, numBendingElements) = node0;
            BendingElements(1, numBendingElements) = node1;
            BendingElements(2, numBendingElements) = node2;
            BendingElements(3, numBendingElements) = node3;
            numBendingElements++;
        }
        else // more than two elements share this edge
        { 
            cout << "Error in finding bending elements." << endl;
            cout << "Found " << indices.size() << " elements that share this edge." << endl;
            return;
        }
    }
    // Matrix4Xi BendingElementsN(4, numBendingElements);
    // for (int i = 0; i < numBendingElements; i++)
    // {
    //     BendingElementsN(0, i) = BendingElements(0, i);
    //     BendingElementsN(1, i) = BendingElements(1, i);
    //     BendingElementsN(2, i) = BendingElements(2, i);
    //     BendingElementsN(3, i) = BendingElements(3, i);
    // }
    cout << "Number of elements: " << numElements << endl;
    cout << "Number of bending elements: " << numBendingElements << endl;
}

void Simulation::getForceJacobian()
{
    if (currentTime <= totalTime / 2) 
    {
        h_running = h_min + (h_max - h_min) * currentTime / (totalTime / 2);
    }
    else 
    {
        h_running = h_max + (h_min - h_max) * (currentTime - totalTime / 2) / (totalTime / 2);
    }
    kb_running = kb * pow(h_running, 3) / pow(h_max, 3);
    ks_running = ks * h_running / h_max;
    W_running = W * h_running / h_max;

    int numDOF = q.size();
    Fs = VectorXd::Zero(numDOF);
    Js = MatrixXd::Zero(numDOF, numDOF);
    for (int i = 0; i < numEdges; i++)
    {
        int n0 = Edges(0, i);
        int n1 = Edges(1, i);
        Vector3d node0 = {q(3 * n0), q(3 * n0 + 1), q(3 * n0 + 2)};
        Vector3d node1 = {q(3 * n1), q(3 * n1 + 1), q(3 * n1 + 2)};
        gradEs_hessEs_Shell gradEs_hessEs_Shell(node0, node1, refLen(i), ks_running(i));
        VectorXd dFs = gradEs_hessEs_Shell.getdF();
        MatrixXd dJs = gradEs_hessEs_Shell.getdJ();
        VectorXi ind(6);
        ind << 3 * n0, 3 * n0 + 1, 3 * n0 + 2, 3 * n1, 3 * n1 + 1, 3 * n1 + 2;
        Fs(ind) = Fs(ind) + dFs;
        Js(ind, ind) = Js(ind, ind) + dJs;
    }
    Fb = VectorXd::Zero(numDOF);
    Jb = MatrixXd::Zero(numDOF, numDOF);
    for (int i = 0; i < numEdges; i++) {
        int n0 = BendingElements(0, i);
        int n1 = BendingElements(1, i);
        int n2 = BendingElements(2, i);
        int n3 = BendingElements(3, i);
        // cout << 3 * n0 << " " << 3 * n1 << " " << 3 * n2 << " " << 3 * n3 << endl;
        Vector3d node0 = {q(3 * n0), q(3 * n0 + 1), q(3 * n0 + 2)};
        Vector3d node1 = {q(3 * n1), q(3 * n1 + 1), q(3 * n1 + 2)};
        Vector3d node2 = {q(3 * n2), q(3 * n2 + 1), q(3 * n2 + 2)};
        Vector3d node3 = {q(3 * n3), q(3 * n3 + 1), q(3 * n3 + 2)};
        double thetaBar = 0;
        gradEb_hessEb_Shell gradEb_hessEb_Shell(node0, node1, node2, node3, thetaBar, kb_running);
        VectorXd dFb = gradEb_hessEb_Shell.getdF();
        MatrixXd dJb = gradEb_hessEb_Shell.getdJ();
        VectorXi ind(12);
        ind << 3 * n0, 3 * n0 + 1, 3 * n0 + 2, 3 * n1, 3 * n1 + 1, 3 * n1 + 2, 3 * n2, 3 * n2 + 1, 3 * n2 + 2, 3 * n3, 3 * n3 + 1, 3 * n3 + 2;
        Fb(ind) = Fb(ind) + dFb;
        Jb(ind, ind) = Jb(ind, ind) + dJb;
    }

    Fv = VectorXd::Zero(numDOF);
    Jv = MatrixXd::Zero(numDOF, numDOF);
    Fv = -visc * (q - q0) / dt * viscLen;
    Jv = -visc * (viscLen / dt) * MatrixXd::Identity(numDOF, numDOF);

    VectorXd elementwise_product = (MassVector / dt).array() * ((q - q0) / dt - qDot0).array();
    Force = elementwise_product + Fs + Fb - Fv - W_running;
    Jac = (MassMatrix / dt) * (MassMatrix / dt) + Js + Jb - Jv;
}

void Simulation::run()
{
    Vector3d g = {0, 0, -9.81};
    int plotStep = 1;

    generateMesh();

    q = VectorXd(3 * numNodes);
    for (int i = 0; i < numNodes; i++)
    {
        q(3 * i) = Nodes(0, i);
        q(3 * i + 1) = Nodes(1, i);
        q(3 * i + 2) = Nodes(2, i);
    }

    double totalM = l * w * h_max * rho;
    double dm = totalM / numNodes;
    MassVector = VectorXd::Constant(3 * numNodes, dm);
    MassMatrix = MatrixXd::Zero(3 * numNodes, 3 * numNodes);
    for (int i = 0; i < 3 * numNodes; i++) {
        MassMatrix(i, i) = MassVector(i);
    }

    W = VectorXd(3 * numNodes);
    for (int i = 0; i < numNodes; i++)
    {
        W(3 * i) = dm * g(0);
        W(3 * i + 1) = dm * g(1);
        W(3 * i + 2) = dm * g(2);
    }

    refLen = VectorXd::Zero(numEdges);
    viscLen = 0.0;
    for (int i = 0; i < numEdges; i++)
    {
        int n0 = Edges(0, i);
        int n1 = Edges(1, i);
        Vector3d diff;
        diff << Nodes(0, n0) - Nodes(0, n1), Nodes(1, n0) - Nodes(1, n1), Nodes(2, n0) - Nodes(2, n1);
        refLen(i) = diff.norm();
        viscLen += refLen(i);
    }
    viscLen /= numEdges;

    ks = 0.5 * sqrt(3) * Y * h_max * refLen.array().pow(2);
    kb = 2 / sqrt(3) * Y * pow(h_max, 3) / 12;
    double tolerance = kb / l * 0.01;

    VectorXi freeDOF(3 * numNodes - 3);
    for (int i = 0; i < 3 * numNodes - 3; i++) 
    {
        freeDOF(i) = i + 3;
    }
    int numFreeDOF = freeDOF.size();
    qDot = q * 0;
    qDot0 = qDot;

    int numSteps = floor(totalTime / dt);
    int maxIter = 20;

    for (int step = 0; step < numSteps; step++)
    {
        double currentTime = dt * step;
        q0 = q;
        qDot0 = qDot;
        VectorXd qFree = q0(freeDOF);
        double error = 10 * tolerance;
        cout << "Initial error value: " << error << endl;
        int iter = 0;
        while (error > tolerance && iter < maxIter)
        {
            getForceJacobian();
            VectorXd forceFree = Force(freeDOF);
            for (int i = 0; i < numFreeDOF; i++) {
                cout << forceFree(i) << " ";
            }
            MatrixXd jacFree = Jac(freeDOF, freeDOF);

            // TODO: (option 1) sparse matrix mldivide via Eigen

            // Eigen::SparseMatrix<double> A(numFreeDOF, numFreeDOF);
            // for (int i = 0; i < numFreeDOF; i++) {
            //     for (int j = 0; j < numFreeDOF; j++) {
            //         A(i, j) = jacFree(i, j);
            //     }
            // }
            // SparseLU<Eigen::SparseMatrix<double>> solver;
            // solver.compute(A);
            // if (solver.info() != Success) 
            // {
            //     cout << "Eigen::SparseLU: decomposition failed." << endl;
            //     return;
            // }
            // dq = solver.solve(Force);
            // if (solver.info() != Success) 
            // {
            //     cout << "Eigen::SparseLU: solving failed." << endl;
            //     return;
            // }

            // TODO: (option 2) use Pardiso to solve linear system
            // dq = JacFree \ ForceFree; 

            // current version: use Eigen/Dense to solve, but may be incorrect and slow
            dq = jacFree.householderQr().solve(forceFree);
            qFree = qFree - dq;
            for (int i = 0; i < numFreeDOF; i++) 
            {
                q(i + 3) = qFree(i);
            }

            error = forceFree.array().abs().sum();
            iter++;
            cout << "Error: " << error << endl;
        }

        if (iter >= maxIter)
        {
            cout << "Could not converge." << endl;
            return;
        }

        VectorXd qDotFree = (qFree - q0(freeDOF)) / dt;
        qDot(freeDOF) = qDotFree;
        
        // write and save mesh
    }
}