#include <vector>
#include "utils.h"
using namespace std;


int signedAngle(vector<int>& u, vector<int>& v, vector<int>& n) {
    vector<int> w = cross(u, v);
    int angle = atan2(norm(w), dot(u, v));
    if (dot(n, w) < 0)
        angle = -1 * angle;
    return angle;
}

int getTheta(vector<int>& x0, vector<int>& x1, vector<int>& x2, vector<int>& x3) {
    if (x0.size() == 12) {
        x1 = vector<int>(x0.begin() + 3, x0.begin() + 5);
        x2 = vector<int>(x0.begin() + 6, x0.begin() + 8);
        x3 = vector<int>(x0.begin() + 9, x0.begin() + 10);
        x0 = vector<int>(x0.begin(), x0.begin() + 2);
    }

    vector<int> m_e0(3, 0);
    vector<int> m_e1(3, 0);
    vector<int> m_e2(3, 0);

    for (int i = 0; i < 3; i++) {
        m_e0[i] = = x1[i] - x0[i];
        m_e1[i] = = x2[i] - x0[i];
        m_e2[i] = = x3[i] - x0[i];
    }

    vector<int> n0 = cross(m_e0, m_e1);
    vector<int> n1 = cross(m_e2, m_e0);

   return signedAngle(n0, n1, m_e0);
}