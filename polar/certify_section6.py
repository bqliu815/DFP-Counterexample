# SPDX-License-Identifier: MIT
"""Certify the Section 6 finite interpolant from its five stored binary64 columns.

Usage: python3 polar/certify_section6.py path/to/matched_finite.mat [--output FILE]
The default output is certificate.json beside the input. Requires NumPy and h5py.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from math import comb
from pathlib import Path

import h5py
import numpy as np


def verify_cutoff():
    """Certify the uniform perturbation bound (21/2) ||a||/sigma.

    In the support annulus put r=||x-x_k||/sigma and z=(3r-1)/2. Write
    a/||a||=(d,e) in the orthonormal radial and tangent frame. Directly
    differentiating the cutoff times a^T(x-x_k) and factoring ||a||/sigma
    gives the symmetric matrix [[d*A,e*B],[e*B,d*B]], where B=phi'(r)
    and A=2*phi'(r)+r*phi''(r). The two absolute row sums are bounded by
    sqrt(A^2+B^2) and sqrt(2)*|B| by Cauchy--Schwarz. In the constant
    cutoff regions the perturbation vanishes.
    """
    a = [0, -45, -45, 360, -270]  # -45*z*(1-z)*(1+2*z-6*z^2)
    b = [0, 0, -45, 90, -45]       # -45*z^2*(1-z)^2
    q = [sum(poly[i] * poly[s-i] for poly in (a, b)
             for i in range(5) if 0 <= s-i <= 4) for s in range(9)]
    max_coefficients = []
    for panel in range(4):
        # On the panel z=(panel+w)/4, convert A^2+B^2 from powers of w
        # to the Bernstein basis of degree eight. This basis is nonnegative
        # and sums to one on [0,1], so its largest coefficient is an upper bound.
        power = [sum(F(q[n]*comb(n,k)*panel**(n-k), 4**n)
                     for n in range(k,9)) for k in range(9)]
        bernstein = [sum(power[k]*F(comb(i,k),comb(8,k)) for k in range(i+1))
                     for i in range(9)]
        max_coefficients.append(max(bernstein))
    constant = F(21,2)
    assert max(max_coefficients) < constant**2
    assert 2*F(45,16)**2 < constant**2  # z(1-z) <= 1/4
    return constant, max_coefficients


def dyadic_integers(array):
    """Represent every stored value as a single integer times 2**(-scale)."""
    nonzero = np.abs(array[array != 0])
    scale = int(53 - np.frexp(nonzero)[1].min())
    def convert(x):
        n, d = float(x).as_integer_ratio()
        return n << (scale-(d.bit_length()-1))
    return [tuple(map(convert, row)) for row in array], scale


def certify(path):
    if not __debug__:
        raise RuntimeError("The exact certificate requires Python assertions. Run without -O.")
    with h5py.File(path, "r") as f:
        data = np.asarray(f["fixedData"]).T
    assert data.dtype == np.float64 and data.shape == (200001, 5) and np.isfinite(data).all()
    constant, coefficients = verify_cutoff()
    ratio_bound = F(503,6250)    # 0.08048
    rows, exponent = dyadic_integers(data)
    assert min(row[4] for row in rows) > 0
    cell_width = 2*max(row[4] for row in rows)
    # Balls in grid cells separated by at least one complete intervening cell
    # have center distance > cell_width >= sum of their radii. Check the
    # remaining candidate pairs exactly once with integer squared distances.
    buckets = {}
    candidates = 0
    for x, y, ax, ay, radius in rows:
        assert ((ax*ax + ay*ay)*ratio_bound.denominator**2
                <= (ratio_bound.numerator*radius)**2)
        u, v = x//cell_width, y//cell_width
        for i in range(u-1,u+2):
            for j in range(v-1,v+2):
                for px, py, pr in buckets.get((i,j), ()):
                    assert (x-px)**2+(y-py)**2 > (radius+pr)**2
                    candidates += 1
        buckets.setdefault((u,v),[]).append((x,y,radius))
    perturbation = constant*ratio_bound
    assert perturbation < 1
    return dict(status="exact_dyadic_checks_passed",
                data_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                number_of_sites=len(rows), common_dyadic_denominator_exponent=exponent,
                separation="Closed disks checked with exact integer squared distances and an exact spatial grid.",
                candidate_pairs=candidates,
                radial_polynomial_panel_bounds=list(map(str,coefficients)),
                cutoff_hessian_constant=str(constant),
                max_correction_ratio_upper=str(ratio_bound),
                perturbation_upper=str(perturbation),
                hessian_interval_exact=list(map(str,(1-perturbation,1+perturbation))),
                hessian_interval=list(map(float,(1-perturbation,1+perturbation))))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matfile", type=Path)
    parser.add_argument("--output", type=Path, help="certificate JSON destination")
    args = parser.parse_args()
    output = args.output or args.matfile.with_name("certificate.json")
    if output.resolve() == args.matfile.resolve():
        parser.error("Choose an output file different from the input MAT file.")
    result = certify(args.matfile)
    report = json.dumps(result, indent=2)
    output.write_text(report + "\n", encoding="utf-8")
    print(report)
