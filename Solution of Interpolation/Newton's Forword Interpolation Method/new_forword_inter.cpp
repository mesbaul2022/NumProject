#include <bits/stdc++.h>
using namespace std;
int fact(int x)
{
  if (x == 0)
    return 1;
  else
    return x * fact(x - 1);
}
int main()
{
  cout << "Enter the N ";
  int n;
  cin >> n;
  vector<double> x(n), y(n);
  cout << "Enter x1-x2 and y:" << endl;
  double sum = 0;
  for (int i = 0; i < n; i++)
  {
    int x1, x2, y1;
    cin >> x1 >> x2 >> y1;
    x[i] = x2;
    y[i] = sum + y1;
    sum += y1;
  }
  cout << "Enter the range  to find x1-x2 :";
  int init, final;
  cin >> init >> final;

  vector<vector<double>> dif(n, vector<double>(n, 0));
  for (int i = 0; i < n; i++)
  {
    dif[i][0] = y[i];
  }
  for (int j = 1; j < n; j++)
  {
    for (int i = 0; i < n - j; i++)
    {
      dif[i][j] = dif[i + 1][j - 1] - dif[i][j - 1];
    }
  }
  double h = x[1] - x[0];
  double u = (final - x[0]) / h;
  double result = y[0];
  double temp = 1.0;
  for (int i = 1; i < n; i++)
  {
    temp = temp * (u - i + 1);
    result += (temp * dif[0][i]) / fact(i);
  }
  cout << "value of Y at " << init << "-" << final << " : " << result << endl;
  return 0;
}
