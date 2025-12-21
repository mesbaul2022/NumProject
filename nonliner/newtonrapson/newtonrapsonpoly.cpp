#include <bits/stdc++.h>
using namespace std;
int degree;
vector<double> coff;
void print()
{
    int d = degree;
    bool m = true;
    for (int i = 0; i < degree + 1; i++)
    {
        if (m)
        {
            cout << coff[i] << "X^" << d;
            m = false;
            d--;
            continue;
        }
        if (coff[i] == 0)
        {
            d--;
            continue;
        }

        if (coff[i] > 0)
            cout << "+";
        if (i != degree)
            cout << coff[i] << "X^" << d;
        else
            cout << coff[i];
        d--;
    }
    cout << "=0" << endl;
}
double funcvalue(double x)
{

    double val = 0;
    double deg = coff.size() - 1;
    for (int i = 0; i < coff.size(); i++)
    {
        val += coff[i] * pow(x, deg);
        deg--;
    }
    return val;
}
double derivative(double x)
{
    double val = 0;
    double deg = coff.size() - 1;
    for (int i = 0; i < coff.size() - 1; i++)
    {
        val += coff[i] * deg * pow(x, deg - 1);
        deg--;
    }
    return val;
}
int main()
{
    cout << "degree:";
    cin >> degree;
    cout << "cofficient:";
    coff.resize(degree + 1);
    for (int i = 0; i < degree + 1; i++)
    {
        cin >> coff[i];
    }
    cout << "equation  ";
    print();
    double xmax = sqrt((coff[1] / coff[0]) * (coff[1] / coff[0]) - 2 * (coff[2] / coff[0]));
    double c = -xmax;
    double num = 0;
    double e = 0.00001;
    while (c <= xmax)
    {
        double x0 = c, x1 = c + 0.65;
        double f0 = funcvalue(x0), f1 = funcvalue(x1);
        if (f0 * f1 < 0)
        {
            num++;
            cout << " interval is (" << x0 << "," << x1 << ")" << endl;
            int it = 0;
            do
            {
                double derivative1 = derivative(x1);
                double x2 = x1 - (f1 / derivative1);
                double fx2 = funcvalue(x2);
                it++;
                if (abs(x2 - x1) < e || abs(fx2 - f1) < e)
                {
                    cout << "root " << num << " is: " << x2 << endl;
                    cout << "iteration number:" << it << endl
                         << endl;
                    break;
                }
                x1 = x2;
                f1 = fx2;

            } while (1);
        }
        c = c + 0.65;
    }

    return 0;
}
