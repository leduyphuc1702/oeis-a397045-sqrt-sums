# Counting Distinct Sums of Square Roots: Proof of the Asymptotic Conjecture in OEIS A397045

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24521851.svg)](https://doi.org/10.5281/zenodo.24521851)

**Preprint (Zenodo DOI):** [https://doi.org/10.5281/zenodo.24521851](https://doi.org/10.5281/zenodo.24521851)

Let $a(n)$ be the number of distinct real numbers in $[n,n+1)$ that are finite sums of square roots of positive integers ([A397045](https://oeis.org/A397045)). Charles R. Greathouse IV conjectured $\log a(n)\sim\kappa n^{2/3}$, possibly with $\kappa=(81\zeta(3)/\pi^2)^{1/3}$. The manuscript proves

$$\log a(n)=\Bigl(\frac{81\zeta(3)}{\pi^2}\Bigr)^{1/3}n^{2/3}+O(n^{1/3}),\qquad \Bigl(\frac{81\zeta(3)}{\pi^2}\Bigr)^{1/3}=2.1447\ldots$$

The key step identifies $a(n)$ with a count of partitions into parts $\sqrt d$, $d>1$ squarefree, using the linear independence of square roots of distinct squarefree integers.

## Contents

- `paper.tex`, `paper.pdf`: the manuscript.
- `verify_a397045.py`: exact integer enumeration of $a(1),\dots,a(20)$ against the OEIS data, and a saddle-point comparison up to $n=35$.

## Reproducing

```bash
pip install numpy mpmath
python3 verify_a397045.py      # prints ALL CHECKS PASSED
```

## Status

Preprint, not yet peer-reviewed. Parts of the derivations and code were developed with AI assistance; every computational claim is backed by the script above.

## License

Manuscript text: CC BY 4.0. Code: MIT License.
