#include <bits/stdc++.h>
using namespace std;

const float EPS = 1e-6;

int main() {
    int t;
    cout << "Enter number of test cases: ";
    cin >> t;

    for (int i = 1; i <= t; i++) {
        int n;
        cout << "Enter number of variables " << i << ": ";
        cin >> n;

        vector<vector<float>> matrix(n, vector<float>(n + 1));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j <= n; j++) {
                cin >> matrix[i][j];
            }
        }

        for (int i = 0; i < n; i++) {
            int pivotrow = i;
            while (pivotrow < n && abs(matrix[pivotrow][i]) < EPS) {
                pivotrow++;
            }

            if (pivotrow == n) continue;
            swap(matrix[i], matrix[pivotrow]);

            float pivot = matrix[i][i];
            for (int j = i; j <= n; j++) {
                matrix[i][j] /= pivot;
            }

            for (int k = 0; k < n; k++) {
                if (k != i) {
                    float factor = matrix[k][i];
                    for (int j = i; j <= n; j++) {
                        matrix[k][j] -= factor * matrix[i][j];
                    }
                }
            }
        }
        bool nosolution = false;
        bool infinitesolution = false;
        int rank = 0;

        for (int i = 0; i < n; i++) {
            bool allzero = true;
            for (int j = 0; j < n; j++) {
                if (abs(matrix[i][j]) > EPS) {
                    allzero = false;
                    break;
                }
            }

            if (allzero) {
                if (abs(matrix[i][n]) > EPS) {
                    nosolution = true;
                    break;
                } else {
                    infinitesolution = true;
                }
            } else {
                rank++;
            }
        }

        cout << "Case " << i << ":" << endl;
        if (nosolution) {
            cout << "No Solution" << endl;
        } else if (infinitesolution || rank < n) {
            cout << "Infinite Solutions" << endl;
        } else {
            cout << "Unique Solution:" << endl;
            for (int i = 0; i < n; i++) {
                cout<< matrix[i][n] <<" ";
            }
        }
        cout << endl;

    }

    return 0;
}
