#include<bits/stdc++.h>
using namespace std;

long long fact(int n)
{
    if(n == 0 || n == 1) return 1;
    return n * fact(n - 1);
}

double g(double x)
{
    return exp(x) + x*x*x + sin(x);
}

double g_prime(double x)
{
    return exp(x) + 3*x*x + cos(x);
}

double g_double_prime(double x)
{
    return exp(x) + 6*x - sin(x);
}

vector<vector<double>> buildDiffTable(vector<double>& values)
{
    int sz = values.size();
    vector<vector<double>> table(sz, vector<double>(sz, 0.0));

    for(int i = 0; i < sz; i++)
        table[i][0] = values[i];

    for(int j = 1; j < sz; j++)
        for(int i = 0; i < sz - j; i++)
            table[i][j] = table[i + 1][j - 1] - table[i][j - 1];

    return table;
}

void process(int caseNo, ifstream& fin, ofstream& fout)
{
    int intervalCount;
    fin >> intervalCount;

    double startPoint, endPoint;
    fin >> startPoint >> endPoint;

    double evalPoint;
    fin >> evalPoint;

    double step = (endPoint - startPoint) / intervalCount;

    vector<double> grid(intervalCount), funcValues(intervalCount);
    for(int i = 0; i < intervalCount; i++)
    {
        grid[i] = startPoint + i * step;
        funcValues[i] = g(grid[i]);
    }

    vector<vector<double>> diff = buildDiffTable(funcValues);

    double u = (evalPoint - grid[0]) / step;

    double approxFirst =
        ( diff[0][1]
        + (2*u - 1) * diff[0][2] / fact(2)
        + (3*u*u - 6*u + 2) * diff[0][3] / fact(3)
        ) / step;

    double approxSecond =
        ( diff[0][2]
        + (u - 1) * diff[0][3]
        ) / (step * step);

    double exactFirst = g_prime(evalPoint);
    double exactSecond = g_double_prime(evalPoint);

    double errorFirst = fabs((exactFirst - approxFirst) / exactFirst) * 100.0;
    double errorSecond = fabs((exactSecond - approxSecond) / exactSecond) * 100.0;

    cout << "\nTEST CASE #" << caseNo << "\n";
    cout << "Difference Table:\n";
    for(int i = 0; i < intervalCount; i++)
    {
        for(int j = 0; j < intervalCount; j++)
            cout << setw(12) << diff[i][j];
        cout << "\n";
    }

    cout << fixed << setprecision(6);
    cout << "Numerical g'(x)  = " << approxFirst << "\n";
    cout << "Exact g'(x)      = " << exactFirst << "\n";
    cout << "Numerical g''(x) = " << approxSecond << "\n";
    cout << "Exact g''(x)     = " << exactSecond << "\n";
    cout << "Error in g'(x)   = " << errorFirst << "%\n";
    cout << "Error in g''(x)  = " << errorSecond << "%\n";

    fout << "\nTEST CASE #" << caseNo << "\n";
    fout << fixed << setprecision(6);
    fout << "Numerical g'(x)  = " << approxFirst << "\n";
    fout << "Exact g'(x)      = " << exactFirst << "\n";
    fout << "Numerical g''(x) = " << approxSecond << "\n";
    fout << "Exact g''(x)     = " << exactSecond << "\n";
    fout << "Error in g'(x)   = " << errorFirst << "%\n";
    fout << "Error in g''(x)  = " << errorSecond << "%\n";
}

int main()
{


    ifstream fin("input.txt");
    ofstream fout("output.txt");



    int totalCases;
    fin >> totalCases;

    cout << "Total Test Cases: " << totalCases << "\n";
    fout << "Total Test Cases: " << totalCases << "\n";

    for(int i = 1; i <= totalCases; i++)
        process(i, fin, fout);

    fin.close();
    fout.close();


    return 0;
}
