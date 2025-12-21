#include<bits/stdc++.h>
using namespace std;

long long fact(int n)
{
    if(n == 0 || n == 1) return 1;
    return n * fact(n - 1);
}

double f(double x)
{
    return x*x + sin(x);
}

double f1(double x)
{
    return 2*x + cos(x);
}

double f2(double x)
{
    return 2 - sin(x);
}

vector<vector<double>> backwardDiffTable(vector<double>& values)
{
    int n = values.size();
    vector<vector<double>> table(n, vector<double>(n, 0.0));

    for(int i = 0; i < n; i++)
        table[i][0] = values[i];

    for(int j = 1; j < n; j++)
        for(int i = n - 1; i >= j; i--)
            table[i][j] = table[i][j - 1] - table[i - 1][j - 1];

    return table;
}

void solve(int tc, ifstream& fin, ofstream& fout)
{
    int n;
    fin >> n;

    double a, b;
    fin >> a >> b;

    double X;
    fin >> X;

    double h = (b - a) / n;

    vector<double> x(n), y(n);
    for(int i = 0; i < n; i++)
    {
        x[i] = a + i * h;
        y[i] = f(x[i]);
    }

    vector<vector<double>> diff = backwardDiffTable(y);

    double u = (X - x[n - 1]) / h;

    double dydx =
        ( diff[n - 1][1]
        + (2*u + 1) * diff[n - 1][2] / fact(2)
        + (3*u*u + 6*u + 2) * diff[n - 1][3] / fact(3)
        ) / h;

    double d2ydx2 =
        ( diff[n - 1][2]
        + (u + 1) * diff[n - 1][3]
        ) / (h * h);

    double exact1 = f1(X);
    double exact2 = f2(X);

    double error1 = fabs((exact1 - dydx) / exact1) * 100.0;
    double error2 = fabs((exact2 - d2ydx2) / exact2) * 100.0;

    cout << "\nTEST CASE #" << tc << "\n";
    cout << "Backward Difference Table:\n";
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
            cout << setw(12) << diff[i][j];
        cout << "\n";
    }

    cout << fixed << setprecision(6);
    cout << "Numerical f'(x)  = " << dydx << "\n";
    cout << "Exact f'(x)      = " << exact1 << "\n";
    cout << "Numerical f''(x) = " << d2ydx2 << "\n";
    cout << "Exact f''(x)     = " << exact2 << "\n";
    cout << "Error in f'(x)   = " << error1 << "%\n";
    cout << "Error in f''(x)  = " << error2 << "%\n";

    fout << "\nTEST CASE #" << tc << "\n";
    fout << fixed << setprecision(6);
    fout << "Numerical f'(x)  = " << dydx << "\n";
    fout << "Exact f'(x)      = " << exact1 << "\n";
    fout << "Numerical f''(x) = " << d2ydx2 << "\n";
    fout << "Exact f''(x)     = " << exact2 << "\n";
    fout << "Error in f'(x)   = " << error1 << "%\n";
    fout << "Error in f''(x)  = " << error2 << "%\n";
}

int main()
{


    ifstream fin("input.txt");
    ofstream fout("output.txt");



    int t;
    fin >> t;

    cout << "Total Test Cases: " << t << "\n";
    fout << "Total Test Cases: " << t << "\n";

    for(int i = 1; i <= t; i++)
        solve(i, fin, fout);

    fin.close();
    fout.close();


    return 0;
}
