"""Exact independent local check of the circular boundary conditions at lambda=41491/100000."""
from math import comb
import numpy as np

PRIME = 10007


def inverse(series):
    result = [pow(int(series[0]), -1, PRIME)]
    for j in range(1, len(series)):
        result.append(-result[0] * sum(series[k] * result[j-k] for k in range(1,j+1)) % PRIME)
    return result


def times(left, right):
    return [sum(left[k] * right[j-k] for k in range(j+1)) % PRIME for j in range(len(left))]


def modular_rank(rows):
    matrix = np.array(rows, dtype=np.int64) % PRIME
    rank = 0
    for col in range(matrix.shape[1]):
        candidates = np.flatnonzero(matrix[rank:,col])
        if not len(candidates):
            continue
        pivot = rank + candidates[0]
        matrix[[rank,pivot]] = matrix[[pivot,rank]]
        row = matrix[rank] * pow(int(matrix[rank,col]), -1, PRIME) % PRIME
        factors = matrix[rank+1:,col].copy()
        matrix[rank+1:] = (matrix[rank+1:] - factors[:,None] * row) % PRIME
        rank += 1
    return rank


def check(degree, parity):
    powers = np.arange(parity, degree+1, 2)
    count = len(powers)
    length = 2*count + 6
    a = 3 * pow(10, -1, PRIME) % PRIME
    lam = 41491 * pow(100000, -1, PRIME) % PRIME
    mu = lam**2 % PRIME
    t = [(-1)**j * comb(j+3,3) % PRIME for j in range(length)]
    field = t.copy()
    field[0] = (field[0] - mu) % PRIME
    inverse_field = inverse(field)
    d = [(-4*mu*q) % PRIME for q in times(t, times(inverse_field, inverse_field))]
    d[0] = (d[0]+1) % PRIME
    identity = np.eye(2*count, dtype=np.int64)
    pressure = [identity[:count]]
    velocity = [identity[count:]]
    for j in range(length-1):
        pnext = sum(d[k]*velocity[j-k] for k in range(j+1))
        unext = sum((-1)**k * ((k+1)*pressure[j-k] - velocity[j-k]) for k in range(j+1))
        unext[:-1] -= (powers[1:] * (powers[1:]-1))[:,None] * pressure[j][1:]
        factor = pow(j+1, -1, PRIME)
        pressure.append((pnext % PRIME) * factor % PRIME)
        velocity.append((unext % PRIME) * factor % PRIME)
    # At the two walls reduce x U + Z F_Z modulo Z² = a² − x², x = R − 1.
    rows = np.zeros((length,2*count), dtype=np.int64)
    for index, power in enumerate(map(int,powers)):
        half = power//2
        for k in range(half+1):
            coeff = (-1)**k * comb(half,k) * pow(a,2*(half-k),PRIME) % PRIME
            for j in range(length-2*k):
                rows[j+2*k] += coeff*power*pressure[j][index]
                if j+2*k+1 < length:
                    rows[j+2*k+1] += coeff*velocity[j][index]
    rank = modular_rank(rows)
    assert rank == 2*count, (degree,parity,rank,2*count)
    return rank


if __name__ == '__main__':
    for degree in (4,8,16,32,64):
        print(degree, check(degree,0), check(degree,1))
