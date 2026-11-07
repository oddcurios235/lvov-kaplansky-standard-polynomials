"""Exact subset-DP for homogeneous cyclic matrix slices."""
import sympy as s


def subset_standard(weights, grades, modulus=None):
    m = len(weights)
    n = len(weights[0])
    dp = [[0] * n for _ in range(1 << m)]
    gs = [0] * (1 << m)
    dp[0] = [1] * n

    for mask in range(1, 1 << m):
        ids = [j for j in range(m) if mask >> j & 1]
        r = len(ids)
        gs[mask] = sum(grades[j] for j in ids) % n
        out = dp[mask]
        for pos, j in enumerate(ids):
            old = mask ^ (1 << j)
            shift = gs[old]
            sign = (-1) ** (r - 1 - pos)
            for i in range(n):
                out[i] += sign * dp[old][i] * weights[j][(i + shift) % n]
        if modulus is not None:
            dp[mask] = [v % modulus for v in out]
    return dp, gs


def slice_grade(weights, grades, k, modulus=None):
    m = len(weights)
    n = len(weights[0])
    dp, gs = subset_standard(weights, grades, modulus)
    full = (1 << m) - 1
    out = [[0] * n for _ in range(n)]
    for mask in range(1 << m):
        comp = full ^ mask
        inv = sum(
            1
            for a in range(m)
            if mask >> a & 1
            for b in range(m)
            if comp >> b & 1 and a > b
        )
        sign = (-1) ** (m - mask.bit_count() + inv)
        sh = gs[mask]
        for i in range(n):
            out[i][(i + sh) % n] += (
                sign * dp[mask][i] * dp[comp][(i + sh + k) % n]
            )
    if modulus is not None:
        out = [[v % modulus for v in row] for row in out]
    return s.Matrix(out)


def det_mod(mat, p):
    A = [[int(v) % p for v in row] for row in mat.tolist()]
    n = len(A)
    ans = 1
    for j in range(n):
        k = next((k for k in range(j, n) if A[k][j]), None)
        if k is None:
            return 0
        if k != j:
            A[k], A[j] = A[j], A[k]
            ans = -ans
        ans = ans * A[j][j] % p
        inv = pow(A[j][j], -1, p)
        for i in range(j + 1, n):
            c = A[i][j] * inv % p
            if c:
                for k in range(j + 1, n):
                    A[i][k] = (A[i][k] - c * A[j][k]) % p
                A[i][j] = 0
    return ans % p
