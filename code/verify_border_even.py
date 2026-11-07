#!/usr/bin/env python3
"""Construct and check exact finite even-standard border certificates.

Alternating products use a subset recurrence, with direct permutation sums
for selected small cases. Recorded full-size integer minors certify injectivity.
"""
from __future__ import annotations

import argparse
import json
import random
import time
from itertools import permutations
from pathlib import Path

import numpy as np
import sympy as sp


def req(x, msg):
    if not x:
        raise RuntimeError(msg)


def all_standards(A):
    n = A[0].shape[0]
    k = len(A)
    D = [None] * (1 << k)
    D[0] = np.eye(n, dtype=object)
    for mask in range(1, 1 << k):
        V = np.zeros((n, n), dtype=object)
        for i in range(k):
            if mask >> i & 1:
                sg = (-1) ** ((mask >> (i + 1)).bit_count())
                V += sg * (D[mask ^ (1 << i)] @ A[i])
        D[mask] = V
    return D


def direct_std(A):
    n = A[0].shape[0]
    k = len(A)
    res = np.zeros((n, n), dtype=object)
    for w in permutations(range(k)):
        sg = (-1) ** sum(w[a] > w[b] for a in range(k) for b in range(a + 1, k))
        P = np.eye(n, dtype=object)
        for i in w:
            P = P @ A[i]
        res += sg * P
    return res


def invsign(S, T, i):
    ss = [j for j in range(max(i + 1, (S | T).bit_length())) if S >> j & 1]
    tt = [j for j in range(max(i + 1, (S | T).bit_length())) if T >> j & 1]
    v = (
        sum(s > t for s in ss for t in tt)
        + sum(s > i for s in ss)
        + sum(t < i for t in tt)
    )
    return (-1) ** v


def normal_matrix(B, u):
    k = len(B)
    n = B[0].shape[0]
    D = all_standards(B)
    full = (1 << k) - 1
    L = np.zeros((n * n, n * n), dtype=object)
    for S in range(1 << k):
        T = full ^ S
        L += invsign(S, T, k) * np.kron(D[S], D[T].T)
    Cs = []
    for i in range(k):
        M = np.zeros((n, n * n), dtype=object)
        rest = full ^ (1 << i)
        S = rest
        while True:
            T = rest ^ S
            v = D[T] @ u
            M += invsign(S, T, i) * np.kron(D[S], v.reshape(1, n))
            if not S:
                break
            S = (S - 1) & rest
        Cs.append(M)
    return np.vstack([L] + Cs), D[full], L, Cs


def independent_rows(M, prime=1009):
    # Modular arithmetic selects rows; the selected integer minor is then computed exactly.
    b = {}
    chosen = []
    for j, row in enumerate(M.tolist()):
        v = [int(x) % prime for x in row]
        for p, w in sorted(b.items()):
            c = v[p]
            if c:
                v = [(x - c * y) % prime for x, y in zip(v, w)]
        p = next((i for i, x in enumerate(v) if x), None)
        if p is None:
            continue
        c = pow(v[p], -1, prime)
        b[p] = [(x * c) % prime for x in v]
        chosen.append(j)
        if len(chosen) == M.shape[1]:
            break
    return chosen


def matrix(a):
    return [[int(x) for x in row] for row in a.tolist()]


def run_case(r, N, seed):
    k = 2 * r - 1
    rng = random.Random(seed)
    for attempt in range(8):
        B = [
            np.array(
                [[rng.randrange(-2, 3) for _ in range(N)] for __ in range(N)],
                dtype=object,
            )
            for j in range(k)
        ]
        u = np.array([rng.randrange(1, 4) for _ in range(N)], dtype=object)
        M, S, L, Cs = normal_matrix(B, u)
        sd = int(sp.Matrix(S.tolist()).det(method='domain-ge'))
        if sd == 0:
            continue
        rows = independent_rows(M)
        if len(rows) != N * N:
            continue
        MM = sp.Matrix(M[rows, :].tolist())
        det = int(MM.det(method='domain-ge'))
        req(det != 0, 'characteristic-zero maximal minor zero')

        # Verify normal columns against replaced-input evaluations.
        for a in range(N):
            for b in range(N):
                E = np.zeros((N, N), dtype=object)
                E[a, b] = 1
                out = all_standards(B + [E])[-1]
                req(list(out.flat) == list(L[:, a * N + b]), f'L column {a,b}')
                for i in range(k):
                    B1 = B.copy()
                    B1[i] = E
                    val = all_standards(B1)[-1] @ u
                    req(list(val) == list(Cs[i][:, a * N + b]), f'C column {i,a,b}')

        # Compare with direct word enumeration in the smallest cases.
        if k <= 7:
            req(np.array_equal(direct_std(B), S), 'ordinary-word endpoint')
            E = np.zeros((N, N), dtype=object)
            E[0, 1] = 1
            B1 = B.copy()
            B1[0] = E
            req(np.array_equal(direct_std(B1) @ u, Cs[0][:, 1]), 'ordinary-word C0')

        # Check signs of the full bordered evaluation and selected variations.
        n = N + 1
        A = []
        for b in B:
            a = np.zeros((n, n), dtype=object)
            a[:N, :N] = b
            A.append(a)
        Y = np.zeros((n, n), dtype=object)
        Y[:N, N] = u
        fullvalue = all_standards(A + [Y])[-1]
        expected = np.zeros((n, n), dtype=object)
        expected[:N, N] = S @ u
        req(np.array_equal(fullvalue, expected), 'border value')
        for b in range(N):
            V = np.zeros((n, n), dtype=object)
            V[N, b] = 1
            Av = A.copy()
            Av[0] = V
            dv = all_standards(Av + [Y])[-1]
            # tr(P dv) = sum(P_ab dv_ba) = -row_b C0(P) u.
            req(list(dv[:N, :N].T.flat) == list(-Cs[0][b, :]), 'lower normal sign')

        return {
            'r': r, 'degree': 2 * r, 'N': N, 'n': N + 1,
            'seed': seed, 'attempt': attempt,
            'B': [matrix(b) for b in B], 'u': u.tolist(),
            'endpoint_det': str(sd), 'normal_rank': N * N,
            'minor_rows': rows, 'normal_minor': str(det),
            'value': matrix(expected), 'independent_permutation_check': k <= 7,
            'all_normal_columns_reconstructed': True,
            'full_border_signs_checked': True,
        }
    raise RuntimeError(f'No certificate in bounded search r={r},N={N}')


def main():
    pa = argparse.ArgumentParser()
    pa.add_argument('--output', default='border_even.json')
    pa.add_argument('--case', nargs=2, type=int)
    a = pa.parse_args()
    cases = [(4, 5), (4, 6), (4, 7), (5, 6), (5, 7), (5, 8), (5, 9)]
    if a.case:
        selected = tuple(a.case)
        if selected not in cases:
            pa.error('only the seven documented certificate cases are supported')
        cases = [selected]

    rec = {'status': 'running', 'cases': []}
    t = time.monotonic()
    try:
        for r, N in cases:
            t1 = time.monotonic()
            x = run_case(r, N, 712000 + 100 * r + N)
            rec['cases'].append(x)
            Path(a.output).parent.mkdir(exist_ok=True, parents=True)
            Path(a.output).write_text(json.dumps(rec, indent=2))
            print('s', 2 * r, 'M', N + 1, 'rank', x['normal_rank'],
                  'seconds', time.monotonic() - t1, flush=True)
        rec['status'] = 'passed'
    except Exception as e:
        rec.update(status='failed', error=repr(e))
        raise
    finally:
        rec['seconds'] = time.monotonic() - t
        Path(a.output).parent.mkdir(exist_ok=True, parents=True)
        Path(a.output).write_text(json.dumps(rec, indent=2))
    print('PASSED', rec['seconds'], flush=True)


if __name__ == '__main__':
    main()
