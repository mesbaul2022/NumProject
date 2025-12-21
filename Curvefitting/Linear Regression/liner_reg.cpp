#include <bits/stdc++.h>
using namespace std;
// Linear Regression: y = a + b*x
int main()
{
    int n;
    cout << "Enter the number of data: ";
    cin >> n;

    vector<double> x(n), y(n);

    cout << "Enter the values of x: ";
    for (int i = 0; i < n; i++)
    {
        cin >> x[i];
    }

    cout << "Enter the values of y: ";
    for (int i = 0; i < n; i++)
    {
        cin >> y[i];
    }

    double sx = 0, sy = 0, sxy = 0, sx_2 = 0;

    for (int i = 0; i < n; i++)
    {
        sx += x[i];
        sy += y[i];
        sxy += x[i] * y[i];
        sx_2 += x[i] * x[i];
    }
    double b = (n * sxy - sx * sy) / (n * sx_2 - sx * sx);
    double a = (sy - b * sx) / n;
    cout << fixed << setprecision(4);
    cout << "The linear equation is: y = " << a << " + " << b << "*x" << endl;
    double x_val = 5;
    cout << "The value of y at x = " << x_val << " is: " << (a + b * x_val) << endl;
    return 0;
}
