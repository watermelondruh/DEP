#include <bits/stdc++.h>
#include "utils.h"
using namespace std;

vector<float> cross(vector<float>& A, vector<float>& B) {
    vector<float> C(3, 0)
    C[0] = A[1] * B[2] - A[2] * B[1];
    C[1] = A[2] * B[0] - A[0] * B[2];
    C[2] = A[0] * B[1] - A[1] * B[0];
    return C;
}

float dot(vector<float>& A, vector<float>& B) {
    return A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
}

float norm(vector<float> A) {
    return pow(A[0] * A[0] + A[1] * A[1] + A[2] * A[2], 0.5);
}

float mean(vector<float> A) {
    int n = A.size();
    float sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += A[i];
    }
    sum /= n;
    return sum;
}

vector<vector<float>> transpose(<vector<vector<float>>& V) {
    float rows = V.size(), cols = V[0].size();
    vector<vector<float>> VV(cols, vector<float>(rows));        
    for (int i = 0; i < rows; ++ i) {
        for (int j = 0; j < cols; ++ j) {
            VV[j][i] = V[i][j];
        }
    }        
    return VV;
}

vector<vector<float>> zeros(int sz1, int sz2) {
    vector<vector<float>> Z(sz1, vector<float>(sz2, 0.));
    return Z;
}

vector<vector<float>> zeros(int sz) {
    vector<vector<float>> Z(sz, vector<float>(sz, 0.));
    return Z;
}




