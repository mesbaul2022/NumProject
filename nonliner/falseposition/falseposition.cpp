#include <bits/stdc++.h>

using namespace std;
double it = 0;
double func(double x)
{
    return pow(x, 4) - 5 * pow(x, 2) + 4;
}
double Xmax(vector<double> &cof, int n)
{
    // return sqrt(pow(an_1 / an, 2) - 2 * (an_2 / an));
    return sqrt(pow(cof[1] / cof[0], 2) - 2 * (cof[2] / cof[0]));
}
double mid(double x1, double x2)
{
    return (x1 * func(x2) - x2 * func(x1)) / (func(x2) - func(x1));
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
        double x1 = xmin;
        double x2 = i;
        double f1 = func(x1);
        double f2 = func(x2);
        if (f1 * f2 < 0)
        {
            it = 0;
            while (true)
            {
                it++;
                x = mid(x1, x2);
                if (func(x1) * func(x) < 0)
                {
                    x2 = x;
                }
                else
                {
                    x1 = x;
                }
                if (fabs(func(x)) < 1e-4 || fabs(x1 - x2) < 1e-4)
                {
                    cout << "Value :" << x << endl;
                    cout << "Number of itretion :" << it << endl;
                    break;
                }
            }
            xmin = x + step;
            i = xmin + step;
        }
    }
}