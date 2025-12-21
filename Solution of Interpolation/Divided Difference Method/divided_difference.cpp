#include <bits/stdc++.h>
using namespace std;

int main()
{
    int n;
    cout << "Enter the number of value : ";
    cin >> n;
    cout << "Enter the value of x and y :" << endl;
    vector<double> x(n), y(n);
    for (int i = 0; i < n; i++)
    {
        cin >> x[i] >> y[i];
    }
    double value;
    cout << "Enter the value of x : ";
    cin >> value;

    vector<vector<double>> table(n, vector<double>(n, 0));

    for (int i = 0; i < n; i++)
    {
         table[i][0] = y[i];
    }
    
    for (int j = 1; j < n; j++)
    {
        for (int i = 0; i < n - j; i++)
        {
            table[i][j] = (table[i + 1][j - 1] - table[i][j - 1]) / (x[i + j] - x[i]);
        }
    }

    cout << "Divided Difference Table :" << endl;
    for (int j = 0; j < n; j++)
    {
        for (int i = 0; i < n - j; i++)
        {
            cout << setw(10) << table[i][j] << " ";
        }
        cout << endl;
    }

    double res = y[0], term = 1.0;
    for (int i = 1; i < n; i++)
    {
        term *= (value - x[i - 1]);
        res += term * table[0][i];
    }

    cout << "Result = " << res << endl;
    return 0;
}
