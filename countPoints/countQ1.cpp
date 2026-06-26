#include <bits/stdc++.h>
#include <cmath>
using namespace std;

// Computes x^n mod p

bool isPrime(int n) {
  if ((n == 0) || (n == 1) || (n == -1))
    return 0;

  for (int i = 2; i <= sqrt(n); i++)
    if (n % i == 0)
      return 0;
  return 1;
}
int powerMod(int x, int n, int p) {
  int val = x;
  for (int i = 0; i < n - 1; i++)
    val = (val * x) % p;
  return val;
}

int legendreSym(int u, int p) {
  u = u % p;
  if (u == 0)
    return 0;

  int val = 1;
  int base = u;
  int ex = (p - 1) / 2;

  while (ex > 0) {
    if (ex & 1)
      val = val * base % p;
    base = base * base % p;
    ex >>= 1;
  }
  return (val == 1) ? 1 : -1;
}

vector<int> legendreTable(int p) {
  vector<int> result;
  for (int i = 0; i < p; i++) {
    result.push_back(legendreSym(i, p));
  }
  return result;
}

int D(int x1, int x2, int x3, int x4, int p) {
  int val =
      powerMod(x1, 4, p) * powerMod(x2, 2, p) * powerMod(x3, 2, p) -
      2 * powerMod(x1, 2, p) * powerMod(x2, 4, p) * powerMod(x3, 2, p) +
      powerMod(x2, 6, p) * powerMod(x3, 2, p) +
      4 * powerMod(x1, 2, p) * powerMod(x2, 3, p) * powerMod(x3, 3, p) -
      4 * powerMod(x2, 5, p) * powerMod(x3, 3, p) -
      2 * powerMod(x1, 2, p) * powerMod(x2, 2, p) * powerMod(x3, 4, p) +
      6 * powerMod(x2, 4, p) * powerMod(x3, 4, p) -
      4 * powerMod(x2, 3, p) * powerMod(x3, 5, p) +
      powerMod(x2, 2, p) * powerMod(x3, 6, p) -
      4 * powerMod(x1, 4, p) * powerMod(x2, 2, p) * x3 * x4 +
      6 * powerMod(x1, 3, p) * powerMod(x2, 3, p) * x3 * x4 -
      2 * powerMod(x1, 2, p) * powerMod(x2, 4, p) * x3 * x4 +
      2 * x1 * powerMod(x2, 5, p) * x3 * x4 - 2 * powerMod(x2, 6, p) * x3 * x4 -
      6 * powerMod(x1, 3, p) * powerMod(x2, 2, p) * powerMod(x3, 2, p) * x4 -
      2 * x1 * powerMod(x2, 4, p) * powerMod(x3, 2, p) * x4 +
      8 * powerMod(x2, 5, p) * powerMod(x3, 2, p) * x4 +
      2 * powerMod(x1, 2, p) * powerMod(x2, 2, p) * powerMod(x3, 3, p) * x4 -
      2 * x1 * powerMod(x2, 3, p) * powerMod(x3, 3, p) * x4 -
      12 * powerMod(x2, 4, p) * powerMod(x3, 3, p) * x4 +
      2 * x1 * powerMod(x2, 2, p) * powerMod(x3, 4, p) * x4 +
      8 * powerMod(x2, 3, p) * powerMod(x3, 4, p) * x4 -
      2 * powerMod(x2, 2, p) * powerMod(x3, 5, p) * x4 +
      4 * powerMod(x1, 4, p) * powerMod(x2, 2, p) * powerMod(x4, 2, p) -
      8 * powerMod(x1, 3, p) * powerMod(x2, 3, p) * powerMod(x4, 2, p) +
      5 * powerMod(x1, 2, p) * powerMod(x2, 4, p) * powerMod(x4, 2, p) -
      2 * x1 * powerMod(x2, 5, p) * powerMod(x4, 2, p) +
      powerMod(x2, 6, p) * powerMod(x4, 2, p) +
      2 * powerMod(x1, 4, p) * x2 * x3 * powerMod(x4, 2, p) +
      4 * powerMod(x1, 3, p) * powerMod(x2, 2, p) * x3 * powerMod(x4, 2, p) -
      2 * x1 * powerMod(x2, 4, p) * x3 * powerMod(x4, 2, p) -
      4 * powerMod(x2, 5, p) * x3 * powerMod(x4, 2, p) +
      4 * powerMod(x1, 3, p) * x2 * powerMod(x3, 2, p) * powerMod(x4, 2, p) -
      7 * powerMod(x1, 2, p) * powerMod(x2, 2, p) * powerMod(x3, 2, p) *
          powerMod(x4, 2, p) +
      10 * x1 * powerMod(x2, 3, p) * powerMod(x3, 2, p) * powerMod(x4, 2, p) +
      6 * powerMod(x2, 4, p) * powerMod(x3, 2, p) * powerMod(x4, 2, p) +
      2 * powerMod(x1, 2, p) * x2 * powerMod(x3, 3, p) * powerMod(x4, 2, p) -
      6 * x1 * powerMod(x2, 2, p) * powerMod(x3, 3, p) * powerMod(x4, 2, p) -
      4 * powerMod(x2, 3, p) * powerMod(x3, 3, p) * powerMod(x4, 2, p) +
      powerMod(x2, 2, p) * powerMod(x3, 4, p) * powerMod(x4, 2, p) -
      4 * powerMod(x1, 4, p) * x2 * powerMod(x4, 3, p) +
      6 * powerMod(x1, 3, p) * powerMod(x2, 2, p) * powerMod(x4, 3, p) -
      6 * powerMod(x1, 2, p) * powerMod(x2, 3, p) * powerMod(x4, 3, p) +
      4 * x1 * powerMod(x2, 4, p) * powerMod(x4, 3, p) -
      6 * powerMod(x1, 3, p) * x2 * x3 * powerMod(x4, 3, p) +
      8 * powerMod(x1, 2, p) * powerMod(x2, 2, p) * x3 * powerMod(x4, 3, p) -
      8 * x1 * powerMod(x2, 3, p) * x3 * powerMod(x4, 3, p) -
      2 * powerMod(x1, 2, p) * x2 * powerMod(x3, 2, p) * powerMod(x4, 3, p) +
      4 * x1 * powerMod(x2, 2, p) * powerMod(x3, 2, p) * powerMod(x4, 3, p) +
      powerMod(x1, 4, p) * powerMod(x4, 4, p);
  return val % p;
}

int main(int argc, char *argv[]) {

  if (argc != 2) {
    cerr << "usage: " << argv[0] << " prime\n";
    return 1;
  }

  int upperBound;

  stringstream convert{argv[1]};

  if (!(convert >> upperBound)) {
    cerr << "usage: " << argv[0] << " prime\n";
    return 1;
  }

  for (int p = 2; p < 1000; p++) {

    if (!isPrime(p))
      continue;

    auto powers = legendreTable(p);
    int count = 0;
    int x1 = 1;
    for (int x2 = 0; x2 < p; x2++) {
      for (int x3 = 0; x3 < p; x3++) {
        for (int x4 = 0; x4 < p; x4++) {
          int Dx = D(x1, x2, x3, x4, p);
          Dx = (Dx + p) % p;
          int legre = powers[Dx];
          if (legre == -1)
            continue;
          if (legre == 1)
            count = count + 2;
          if (legre == 0)
            count++;
        }
      }
    }

    x1 = 0;
    int x2 = 1;
    for (int x3 = 0; x3 < p; x3++) {
      for (int x4 = 0; x4 < p; x4++) {
        int Dx = D(x1, x2, x3, x4, p);
        Dx = (Dx + p) % p;
        int legre = powers[Dx];
        if (legre == -1)
          continue;
        if (legre == 1)
          count = count + 2;
        if (legre == 0)
          count++;
      }
    }

    x1 = 0;
    x2 = 0;
    int x3 = 1;
    for (int x4 = 0; x4 < p; x4++) {
      int Dx = D(x1, x2, x3, x4, p);
      Dx = (Dx + p) % p;
      int legre = powers[Dx];
      if (legre == -1)
        continue;
      if (legre == 1)
        count = count + 2;
      if (legre == 0)
        count++;
    }

    int Dx = D(0, 0, 0, 1, p);
    Dx = (Dx + p) % p;
    int legre = powers[Dx];
    if (legre == -1)
      continue;
    if (legre == 1)
      count = count + 2;
    if (legre == 0)
      count++;

    cout << "F_" << p << " points for D: " << count << "\n";
  }
}
