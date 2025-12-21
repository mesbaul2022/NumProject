#include <bits/stdc++.h>
using namespace std;
vector<double> cof;
int deg;
double valueFunc(double x)
{
    int n = deg;
    double value = 0;
    for (int i = 0; i < cof.size(); i++)
    {
        value += cof[i] * pow(x, n);
        n--;
    }
    return value;
}

void FunctionPrint()
{
    int n = deg;
    for (int i = 0; i <= deg; i++)
    {
        if (cof[i] == 0)
        {
            n--;
            continue;
        }
        if (i != 0 && cof[i] > 0)
        {
            cout << "+";
        }
        if (cof[i] < 0)
        {
            cout << "-";
        }
        double c = fabs(cof[i]);
        if (!(c == 1 && n != 0))
        {
            cout << c;
        }
        if (n > 0)
        {
            cout << "X";
        }
        if (n > 1)
        {
            cout << "^" << n;
        }
        n--;
    }
    cout << "=0" << endl;
}
double X_i(double x1, double x2)
{
    return (x1 * valueFunc(x2) - x2 * valueFunc(x1)) / (valueFunc(x2) - valueFunc(x1));
}
int main()
{
    cout << "Enter DEGREE ";
    cin >> deg;
    cof.resize(deg + 1);
    cout << "Enter thr coffitient ";
    for (int i = 0; i < deg + 1; i++)
    {
        cin >> cof[i];
    }
    FunctionPrint();
    double xmax = 0;
    for (int i = 0; i < cof.size(); i++)
    {
        xmax = max(xmax, cof[i] / cof[0]);
    }
    xmax += 1;
    double xmin = (-1) * xmax;
    double step = 0.45;
    int num = 0;
    double x2;
    for (double i = xmin + step; i <= xmax; i += step)
    {
        double x0 = xmin;
        double x1 = i;
        double f1 = valueFunc(x0);
        double f2 = valueFunc(x1);
        if (f1 * f2 < 0)
        {
            cout << "The range for root " << ++num << ": " << x0 << "--->" << x1 << endl;
            while (true)
            {
                double temp = x1;
                x1 = X_i(x0, x1);
                x0 = temp;
                if (fabs(x1 - x0) < 1e-3 && fabs(valueFunc(x1)) < 1e-3)
                {
                    cout << "Root " << num << " :" << x1 << endl;
                    break;
                }
            }
            xmin = x1 + step;
            i = xmin + step;
        }
    }
}