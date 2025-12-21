#include <bits/stdc++.h>
using namespace std;

// Function to integrate
float func(float x) {
    float s = sin(x);
    return pow(s, 5) + 4 * pow(s, 4) + 1;
}

int main() {
    float a, b;
    int n;
    cout << "Enter the values of a, b and n (n must be divisible by 3): ";
    cin >> a >> b >> n;
    float h = (b - a) / n;
    float sum = func(a) + func(b);

    for (int i = 1; i < n; i++) {
        float x = a + i * h;
        if (i % 3 == 0) {
            sum += 2 * func(x);
        } else {
            sum += 3 * func(x);
        }
    }

    float result = (3 * h / 8) * sum;
    cout << fixed << setprecision(6);
    cout << "Integral result = " << result << endl;

    return 0;
}
