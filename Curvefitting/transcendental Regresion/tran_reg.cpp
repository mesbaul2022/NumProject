#include <bits/stdc++.h>
using namespace std;

// fitting a linear line
// y = a + be^(x/4)
// y = a + bx' = a+bf(x)
// x' = f(x) = e^(x/4)

pair<double, double> linear_reg(vector<double> vx, vector<double> vy)
{
    double n = vx.size();
    double x = 0, xy = 0, y = 0, x2 = 0;
    for (int i = 0; i < n; i++)
    {
        x += exp(vx[i] / 4);
        y += vy[i];
        xy += exp(vx[i] / 4) * vy[i];
        x2 += exp(vx[i] / 4) * exp(vx[i] / 4);
    }
    double b = (n * xy - x * y) / (n * x2 - x * x);
    double a = (y - b * x) / n;
    return {b, a};
}

int main()
{
    int n;
    cout << "Enter the no of points: ";
    cin >> n;
    vector<double> vx(n), vy(n);
    cout << "Enter values of x: " << endl;
    for (int i = 0; i < n; i++)
    {
        cin >> vx[i];
        cin >> vy[i];
    }

    auto [b, a] = linear_reg(vx, vy);
    cout << "y = " << a << " + " << b << "e^(x/4)" << endl;
    cout << "Enter any value of x: ";
    double x;
    cin >> x;
    cout << "y(" << x << ") = " << a + b * exp(x / 4) << endl;
}
