#!/usr/bin/env python3
"""Saved-input verifier for the even-standard bordered-core certificates.

Every normal-map column is recomputed with a first-letter recurrence.
The endpoint determinant is computed exactly; selected minors are checked
modulo two primes. Explicit exceptions remain active under python -O.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import sympy as sp


def need(c, msg):
    if not c:
        raise RuntimeError(msg)


def alternate_first(A):
    k = len(A)
    n = A[0].shape[0]
    values = [None] * (1 << k)
    values[0] = np.eye(n, dtype=object)
    for mask in range(1, 1 << k):
        total = np.zeros((n, n), dtype=object)
        position = 0
        for i in range(k):
            if mask & (1 << i):
                total += ((-1) ** position) * (A[i] @ values[mask ^ (1 << i)])
                position += 1
        values[mask] = total
    return values[-1]


def modular_det(M, p):
    A = [[int(z) % p for z in row] for row in M]
    n = len(A)
    d = 1
    for j in range(n):
        pivot = next((i for i in range(j, n) if A[i][j]), None)
        if pivot is None:
            return 0
        if pivot != j:
            A[j], A[pivot] = A[pivot], A[j]
            d = -d
        v = A[j][j]
        d = d * v % p
        inv = pow(v, -1, p)
        for i in range(j + 1, n):
            v = A[i][j] * inv % p
            if v:
                for h in range(j + 1, n):
                    A[i][h] = (A[i][h] - v * A[j][h]) % p
            A[i][j] = 0
    return d % p


def verify(c):
    B = [np.array(x, dtype=object) for x in c['B']]
    u = np.array(c['u'], dtype=object)
    N = c['N']
    k = len(B)
    S = alternate_first(B)
    detS = int(sp.Matrix(S.tolist()).det(method='domain-ge'))
    need(detS != 0 and str(detS) == c['endpoint_det'], 'endpoint determinant')
    columns = []
    for a in range(N):
        for b in range(N):
            P = np.zeros((N, N), dtype=object)
            P[a, b] = 1
            LP = alternate_first(B + [P])
            need(sum(LP[i, i] for i in range(N)) == 0, 'trace of L')
            col = list(LP.flat)
            for i in range(k):
                C = B.copy()
                C[i] = P
                col.extend(alternate_first(C) @ u)
            columns.append(col)
    M = np.array(columns, dtype=object).T
    selected = M[c['minor_rows'], :].tolist()
    residues = {}
    for p in (101, 103):
        d = modular_det(selected, p)
        need(d != 0, 'zero reconstructed determinant residue')
        need(d == int(c['normal_minor']) % p, 'saved integer minor mismatch modulo prime')
        residues[str(p)] = d

    # Independent evaluation on the enlarged tuple.
    n = N + 1
    A = []
    for b in B:
        a = np.zeros((n, n), dtype=object)
        a[:N, :N] = b
        A.append(a)
    Y = np.zeros((n, n), dtype=object)
    Y[:N, N] = u
    value = alternate_first(A + [Y])
    expected = np.array(c['value'], dtype=object)
    need(np.array_equal(value, expected), 'actual bordered value')
    need(np.array_equal(value @ value, np.zeros_like(value)), 'square zero')
    need(sp.Matrix(value.tolist()).rank() == 1, 'rank one')
    return {
        'r': c['r'], 'N': N, 'normal_rank': N * N,
        'residues': residues,
        'method': 'first-letter recurrence; all source-unit columns; exact modular determinants; actual bordered value',
        'status': 'passed',
    }


def main():
    pa = argparse.ArgumentParser()
    pa.add_argument('--input', required=True)
    pa.add_argument('--output', required=True)
    a = pa.parse_args()
    record = {'status': 'running', 'cases': []}
    t = time.monotonic()
    try:
        data = json.loads(Path(a.input).read_text())
        need(data['status'] == 'passed', 'input record not passed')
        for c in data['cases']:
            tt = time.monotonic()
            result = verify(c)
            record['cases'].append(result)
            Path(a.output).parent.mkdir(exist_ok=True, parents=True)
            Path(a.output).write_text(json.dumps(record, indent=2))
            print(
                'verified', result['r'], result['N'], result['residues'],
                'seconds', time.monotonic() - tt, flush=True,
            )
        record['status'] = 'passed'
    except Exception as e:
        record.update(status='failed', error=repr(e))
        raise
    finally:
        record['seconds'] = time.monotonic() - t
        Path(a.output).parent.mkdir(exist_ok=True, parents=True)
        Path(a.output).write_text(json.dumps(record, indent=2))
    print('PASSED', record['seconds'], flush=True)


if __name__ == '__main__':
    main()
