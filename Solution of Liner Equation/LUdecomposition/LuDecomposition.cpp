#include <bits/stdc++.h>
using namespace std;

int n;
vector<vector<double>> M, L, U;
vector<double> B, Y, result;
const double EPS = 1e-9;

void LU()
{
    for (int i = 0; i < n; i++)
    {
        for (int k = i; k < n; k++)
        {
            double sum = 0;
            for (int j = 0; j < i; j++)
            {
                sum += L[i][j] * U[j][k];
            }
            U[i][k] = M[i][k] - sum;
            if (abs(U[i][k]) < EPS)
                U[i][k] = 0;
        }

        for (int k = i; k < n; k++)
        {
            if (i == k)
            {
                L[i][i] = 1;
            }
            else
            {
                double sum = 0;
                for (int j = 0; j < i; j++)
                {
                    sum += L[k][j] * U[j][i];
                }
                if (abs(U[i][i]) > EPS)
                {
                    L[k][i] = (M[k][i] - sum) / U[i][i];
                }
                else
                {
                    L[k][i] = 0;
                }
            }
        }
    }
}

void FS()
{
    for (int i = 0; i < n; i++)
    {
        double sum = 0;
        for (int j = 0; j < i; j++)
        {
            sum += L[i][j] * Y[j];
        }
        Y[i] = B[i] - sum;
        if (abs(Y[i]) < EPS)
            Y[i] = 0;
    }
}

void BS()
{
    for (int i = n - 1; i >= 0; i--)
    {
        double sum = 0;
        for (int j = i + 1; j < n; j++)
        {
            sum += U[i][j] * result[j];
        }
        result[i] = (Y[i] - sum) / U[i][i];
    }
}

int main()
{
    int test;
    cout << "Enter the test case :";
    cin >> test;
    int caseNum = 1;
    while (test--)
    {
        cout << "Enter the number of variables: ";
        if (!(cin >> n))
            break;

        M.assign(n, vector<double>(n, 0));
        B.assign(n, 0);
        Y.assign(n, 0);
        result.assign(n, 0);
        L.assign(n, vector<double>(n, 0));
        U.assign(n, vector<double>(n, 0));

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                cin >> M[i][j];
            }
        }
        for (int j = 0; j < n; j++)
        {
            cin >> B[j];
        }

        LU();
        FS();

        bool singular = false;
        bool noSolution = false;
        for (int i = 0; i < n; i++)
        {
            if (abs(U[i][i]) < EPS)
            {
                singular = true;
                if (abs(Y[i]) > EPS)
                {
                    noSolution = true;
                    break;
                }
            }
        }

        cout << "Case " << caseNum++ << ":" << endl;
        if (noSolution)
        {
            cout << "No Solution" << endl;
        }
        else if (singular)
        {
            cout << "Infinite Solutions" << endl;
        }
        else
        {
            cout << "Unique Solution:" << endl;
            BS();
            for (int j = 0; j < n; j++)
            {
                cout << result[j] << (j == n - 1 ? "" : " ");
            }
            cout << endl;
        }
        cout << endl;
    }
    return 0;
}
