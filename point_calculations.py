from itertools import product
from collections.abc import Iterator
from math import isqrt

def primes_upto(n : int) -> list[int]:
    out = []
    for x in range(2, n + 1):
        prime = True
        for q in range(2, isqrt(x) + 1):
            if x % q == 0:
                prime = False
                break
                
        if prime:
            out.append(x)
    return out

def chi5(p : int):
    if p % 5 == 0:
        return 0
    return 1 if p % 5 in (1, 4) else -1

def p3_points(p : int) -> Iterator[tuple[int, ...]]:
    # Normalized projective points [x0:x1:x3:x4] in P^3(F_p).
    # The first nonzero coordinate is set equal to 1.

    for k in range(4):
        for tail in product(range(p), repeat = 3 - k):
            yield (0,) * k + (1,) + tail


def root_count_quad(a : int, b : int, c : int, p : int) -> int:
    a %= p
    b %= p
    c %= p
    if a == 0:
        if b == 0:
            return p if c == 0 else 0
        return 1

    if p == 2:
        return sum(1 for t in range(p) if (a*t*t + b*t + c) % p == 0)
    
    disc = (b*b - 4*a*c)
    if disc == 0:
        return 1

    return 2 if pow(disc, (p - 1)//2, p) == 1 else 0

def count_X(p : int) -> int:
    total = 1 #The point [0:0:1:0:0]
    for x0, x1, x3, x4 in p3_points(p):
        A = (x0 - x1) * (x1 - x3)
        L = (x0 - x4)
        M = (x0 + x1 - x3 - x4)
        K = (x0 - x3) * (x0 - x3)
        C0 = (x0 - x1) * (x0 * x1 - x0 * x3 - x1 * x4 + x4 * x4)

        # f(y, t) = a t^2 + bt + c with t = x2
       
        a = x1 * K
        b = A*L*M - x1 * C0 - x1*x4*K
        c = x1*x4*C0
        total += root_count_quad(a,b,c,p)
    return total

def eta_coeffs(N):
    #eta(q)^4 eta(5q)^4 = q prod (1 - q^n)^4 (1 - q^(5n))^4
    poly = [0] * (N + 1)
    poly[1] = 1

    def multiply_by_factor(poly, step):
        coeffs = [1,-4,6,-4,1]
        new = [0] * (N + 1)
        for i, ai in enumerate(poly):
            if ai:
                for j, cj in enumerate(coeffs):
                    k = i + j * step
                    if k <= N:
                        new[k] += ai * cj

        return new

    for n in range(1, N + 1):
        poly = multiply_by_factor(poly, n)
        if 5 * n <= N:
            poly = multiply_by_factor(poly, 5*n)
    return poly

N = 127
coeff = eta_coeffs(N)
for p in primes_upto(N):
    NX = count_X(p)
    ap = 1 + p**3 + 12*p**2 - (17 + chi5(p))*p - NX
    print(p, NX, ap, coeff[p], ap == coeff[p])
