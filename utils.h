#ifndef UTILS_H
#define UTILS_H

#include <vector>
using namespace std;

vector<float> cross(vector<float>& A, vector<float>& B);

float dot(vector<float>& A, vector<float>& B);

float norm(vector<float> A);

float mean(vector<float> A);

vector<vector<float>> transpose(vector<vector<float>>& V);

#endif