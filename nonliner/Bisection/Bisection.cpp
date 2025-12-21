#include <bits/stdc++.h>

using namespace std;
double it;
double func(double x)
{
    return pow(x, 4) - 5 * pow(x, 2) + 4;
}
double Xmax(vector<double> &cof, int n)
{
    return sqrt(pow(cof[1] / cof[0], 2) - 2 * (cof[2] / cof[0]));  
}
double mid(double x1, double x2)
{
    return (x1 + x2) / 2;
}
int main()
{
    // X^4+0X^3-5X^2+4=0;
    int n;
    cout << "Enter the DEGREE OF The equation ";
    cin >> n;
    cout << "Enter the coffetient hign->low ";
    vector<double> cof(n + 1);
    for (int i = 0; i < n + 1; i++)
    {
        cin >> cof[i];
    }
    double xmax = Xmax(cof, n);
    double xmin = xmax * (-1);
    double step = 0.5;
    double x;
    for (double i = xmin + step; i <= xmax; i += step)
    {
        double x0 = xmin;
        double x1 = i;
        double f1 = func(x0);
        double f2 = func(x1);
        if (f1 * f2 < 0)
        {
            cout << "Bracket " << x0 << " " << x1 << endl;
            it = 0;
            while (true)
            {
                it++;
                x = mid(x0, x1);
                if (func(x0) * func(x) < 0)
                {
                    x1 = x;
                }
                else
                {
                    x0 = x;
                }
                if (fabs(func(x)) < 1e-6 || fabs(x1 - x0) < 1e-6)
                {
                    cout << "Value :" << x << endl;
                    cout << "Number of itretion :" << it << endl;
                    break;
                }
            }
            xmin = x + step;
            i = xmin + step;
        }
        xmin+=step;
    }
}