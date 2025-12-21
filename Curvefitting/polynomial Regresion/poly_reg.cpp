#include <bits/stdc++.h>
using namespace std;

void gaus(vector<vector<double>> &matrix, vector<double> &res)
{
    int n = matrix.size();
    for (int i = 0; i < n; i++)
    {
        int piv = i;
        for (int j = i + 1; j < n; j++)
        {
            if (fabs(matrix[j][i]) > fabs(matrix[piv][i]))
            {
                piv = j;
            }
        }
        swap(matrix[i], matrix[piv]);

        for (int j = i + 1; j < n; j++)
        {
            double r = matrix[j][i] / matrix[i][i];
            for (int k = i; k <= n; k++)
            {
                matrix[j][k] -= r * matrix[i][k];
            }
        }
    }
    res.assign(n, 0);
    for (int i = n - 1; i >= 0; i--)
    {
        res[i] = matrix[i][n];
        for (int j = i + 1; j < n; j++)
        {
            res[i] -= matrix[i][j] * res[j];
        }
        res[i] /= matrix[i][i];
    }
}

int main()
{
    int d, n;
    cout << "Enter the degree: ";
    cin >> d;
    cout << "Enter the number of data: ";
    cin >> n;

    vector<double> x(n), y(n);
    cout << "Enter the value of x: ";
    for (int i = 0; i < n; i++)
        cin >> x[i];

    cout << "Enter the value of y: ";
    for (int i = 0; i < n; i++)
        cin >> y[i];

    vector<vector<double>> matrix(d + 1, vector<double>(d + 2, 0));
    for (int i = 0; i <= d; i++)
    {
        for (int j = 0; j <= d; j++)
        {
            for (int k = 0; k < n; k++)
            {
                matrix[i][j] += pow(x[k], i + j);
            }
        }
        for (int k = 0; k < n; k++)
        {
            matrix[i][d + 1] += y[k] * pow(x[k], i);
        }
    }

    vector<double> res;
    gaus(matrix, res);

    cout << "Polynomial coefficients :";
    for (auto x : res)
        cout << x << " ";
    cout << endl;

    return 0;
}
