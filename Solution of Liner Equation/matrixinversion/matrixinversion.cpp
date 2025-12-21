#include <bits/stdc++.h>
using namespace std;

const double EPS = 1e-12;

double detg(vector<vector<double>> a) {
    int n = (int)a.size();
    double d = 1.0;
    int s = 1;

    for (int i = 0; i < n; i++) {
        int p = i;
        for (int r = i; r < n; r++)
            if (fabs(a[r][i]) > fabs(a[p][i])) p = r;

        if (fabs(a[p][i]) < EPS) return 0.0;
        if (p != i) { swap(a[p], a[i]); s = -s; }

        double piv = a[i][i];
        d *= piv;

        for (int r = i + 1; r < n; r++) {
            double f = a[r][i] / piv;
            for (int c = i; c < n; c++) a[r][c] -= f * a[i][c];
        }
    }
    return d * s;
}

int rk(vector<vector<double>> a) {
    int n = (int)a.size();
    int m = (int)a[0].size();
    int r = 0;

    for (int c = 0; c < m && r < n; c++) {
        int p = r;
        for (int i = r; i < n; i++)
            if (fabs(a[i][c]) > fabs(a[p][c])) p = i;

        if (fabs(a[p][c]) < EPS) continue;
        swap(a[p], a[r]);

        for (int i = r + 1; i < n; i++) {
            double f = a[i][c] / a[r][c];
            for (int j = c; j < m; j++) a[i][j] -= f * a[r][j];
        }
        r++;
    }
    return r;
}

vector<vector<double>> adj(vector<vector<double>>& a) {
    int n = (int)a.size();
    vector<vector<double>> ad(n, vector<double>(n, 0.0));

    if (n == 1) { ad[0][0] = 1.0; return ad; }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            vector<vector<double>> m(n - 1, vector<double>(n - 1));
            int r2 = 0;
            for (int r = 0; r < n; r++) if (r != i) {
                int c2 = 0;
                for (int c = 0; c < n; c++) if (c != j) {
                    m[r2][c2++] = a[r][c];
                }
                r2++;
            }
            double cof = detg(m);
            if ((i + j) & 1) cof = -cof;
            ad[j][i] = cof; 
        }
    }
    return ad;
}

int main() {
    int t;
    cout<<"Enter the test case ";
    cin >> t;

    for (int cs = 1; cs <= t; cs++) {
        int n;
        cin >> n;

        vector<vector<double>> A(n, vector<double>(n));
        vector<double> B(n);

        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                cin >> A[i][j];

        for (int i = 0; i < n; i++) cin >> B[i];

        cout << "Case " << cs << ":\n";

        double d = detg(A);

        if (fabs(d) > EPS) {
            vector<vector<double>> ad = adj(A);
            vector<double> x(n, 0.0);

            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    double inv = ad[i][j] / d;
                    x[i] += inv * B[j];
                }
            }

            cout << "Unique solution: ";
            cout << fixed << setprecision(3);
            for (int i = 0; i < n; i++) cout << x[i] << (i + 1 == n ? '\n' : ' ');
        } else {
            vector<vector<double>> Ab(n, vector<double>(n + 1));
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) Ab[i][j] = A[i][j];
                Ab[i][n] = B[i];
            }

            int rA = rk(A);
            int rAb = rk(Ab);

            if (rA < rAb) cout << "NO solution\n";
            else cout << "INFINITE Solution\n";
        }

        cout << "\n";
    }
    return 0;
}
