"""Verification for "Counting distinct sums of square roots" (OEIS A397045).

1. Exact enumeration of a(1..20) with integer arithmetic only (no floating-point
   comparisons): every radical sqrt(d) is enclosed in [w, w+1)/Q, and the
   enumeration asserts that the floor of every enumerated energy is unambiguous.
2. Comparison of log a(n), n<=35 (terms 21..35 from the OEIS entry), with
   (a) the leading term (81 zeta(3)/pi^2)^(1/3) n^(2/3) of the theorem and
   (b) the saddle-point value t X + log Z(t) - log(t sqrt(2 pi Var_t)),
       X = n+1, mean_t = X  (Remark in the paper; not part of the proof).

Usage: python3 verify_a397045.py
"""
from itertools import accumulate
from math import isqrt, log, pi, sqrt

import mpmath as mp
import numpy as np

OEIS = [3, 7, 17, 35, 76, 148, 295, 564, 1068, 1965, 3579, 6397, 11318, 19774, 34156, 58444,
        99050, 166319, 277066, 458000, 751610, 1224935, 1983456, 3191878, 5106583,
        8124040, 12856700, 20243194, 31718970, 49469766, 76809925, 118750958, 182838429,
        280394529, 428363342]

# ------------------------------------------------------------ 1. enumeration
N = 20
Q = 10**12
sf = [d for d in range(2, (N + 1)**2) if all(d % (p * p) for p in range(2, isqrt(d) + 1))]
w = [isqrt(d * Q * Q) for d in sf]          # w <= Q sqrt(d) < w + 1, increasing in d
hist = [0] * (N + 2)


def visit(start, lo, depth):
    # true scaled energy lies in [lo, lo + depth)
    j = lo // Q
    assert depth == 0 or (lo + depth - 1) // Q == j, "ambiguous floor"
    hist[j] += 1
    for i in range(start, len(w)):
        new = lo + w[i]
        if new >= (N + 1) * Q:
            break
        visit(i, new, depth + 1)


visit(0, 0, 0)
A = list(accumulate(hist[:N + 1]))          # A[n] = #{radical vectors with E < n+1}
assert A[1:] == OEIS[:N], "enumeration disagrees with OEIS"
print("1. exact enumeration: a(n) = A(n+1) reproduces a(1..%d) of A397045" % N)

# --------------------------------------------------------- 2. asymptotics
z3 = float(mp.zeta(3))
B = 12 * z3 / pi**2
kappa = 3 * (B / 4)**(1 / 3)
assert abs(kappa - float(mp.root(81 * mp.zeta(3) / mp.pi**2, 3))) < 1e-12
D = 4 * 10**6                                # squarefree d < D: ample for t >= 0.05
sieve = np.ones(D + 1, bool)
sieve[:2] = False
for p in range(2, isqrt(D) + 1):
    sieve[p * p::p * p] = False
sieve[1] = False
roots = np.sqrt(np.flatnonzero(sieve).astype(float))


def stats(t):
    e = np.exp(-t * roots)
    return -np.log1p(-e).sum(), (roots * e / (1 - e)).sum(), (roots**2 * e / (1 - e)**2).sum()


print("2.  n   log a(n)   saddle value   difference   leading term   log a(n)/n^(2/3)")
for n in (5, 10, 15, 20, 25, 30, 35):
    X = n + 1
    lo, hi = 0.01, 5.0
    for _ in range(80):
        t = (lo + hi) / 2
        lz, m, v = stats(t)
        lo, hi = (t, hi) if m > X else (lo, t)
    pred = t * X + lz - log(t * sqrt(2 * pi * v))
    la = log(OEIS[n - 1])
    print("   %2d  %9.4f  %12.4f  %10.4f  %12.4f  %10.4f" % (n, la, pred, la - pred, kappa * n**(2 / 3),
                                                       la / n**(2 / 3)))
print("   constant (81 zeta(3)/pi^2)^(1/3) =", mp.nstr(mp.root(81 * mp.zeta(3) / mp.pi**2, 3), 15))
print("ALL CHECKS PASSED")
