#include <bits/stdc++.h>
using namespace std;

void solve() {
    int n;
    cout<<"Enter the number of variable ";
    cin >> n;
    vector<vector<double>> a(n, vector<double>(n + 1));
    double EPS = 1e-12;
    cout<<"Enter the matrix "<<endl;
    for (int i = 0; i < n; i++)
        for (int j = 0; j <= n; j++) cin >> a[i][j];

    for (int i = 0; i < n; i++) {
        int pivotrow = i;
            while (pivotrow < n && abs(a[pivotrow][i]) < EPS) {
                pivotrow++;
            }

            if (pivotrow == n) continue;
            swap(a[i], a[pivotrow]);

        double t = a[i][i];
        for (int j = 0; j <= n; j++) a[i][j] /= t;

        for (int j = 0; j < n; j++) {
            if (i != j) {
                double r = a[j][i];
                for (int k = 0; k <= n; k++) a[j][k] -= r * a[i][k];
            }
        }
    }

    bool inf = false, none = false;
    for (int i = 0; i < n; i++) {
        bool zero = true;
        for (int j = 0; j < n; j++) if (fabs(a[i][j]) > EPS) zero = false;
        if (zero) {
            if (fabs(a[i][n]) > EPS) none = true;
            else inf = true;
        }
    }

    if (none) cout << "NO solution" << endl;
    else if (inf) cout << "INFINITE Solution" << endl;
    else {
        cout << "Unique solution: ";
        for (int i = 0; i < n; i++) cout << a[i][n] << " ";
        cout << endl;
    }
}

int main() {
    int t;
    cout << "Enter number of test cases: ";
    cin >> t;
    while (t--) {
        solve();
        cout<<endl;
    }
    return 0;
}
