# PyTest script
# run with: 'py.test -s -v test_amc.py'

import pytest
import amc

def test_sort():
    m = [[0,0],[0,0]]
    r = [[0,0],[0,0]]
    assert amc.sort(m) == r

    m = [[0,0], [1,2]]
    r = [[2,1], [0,0]]
    assert amc.sort(m) == r

    m = [[1,2,3],
         [0,0,0],
         [3,2,1]]
    r = [[1,3,2],
         [3,1,2],
         [0,0,0]]
    assert amc.sort(m) == r

    m = [[1,1,1,1], #s0
         [0,0,0,0], #s1
         [0,0,0,0], #s2
         [1,2,3,4]] #s3
    r = [[1,1,1,1], #s0
         [1,4,2,3], #s3
         [0,0,0,0], #s1
         [0,0,0,0]] #s2
    assert amc.sort(m) == r


def calculate_b(m):
    # B = (I-Q)^-1 * R
    m = amc.sort(m)
    n = amc.normalize(m)
    (q, r) = amc.decompose(n)
    i = amc.identity(len(q))
    s = amc.subtract(i, q)
    v = amc.getMatrixInverse(s)
    b = amc.multiply(v, r)
    return b

def test_b():
    # test 1
    m = [[0,1,0,0,0,1],  # s0, the initial state, goes to s1 and s5 with equal probability
         [4,0,0,3,2,0],  # s1 can become s0, s3, or s4, but with different probabilities
         [0,0,0,0,0,0],  # s2 is terminal, and unreachable (never observed in practice)
         [0,0,0,0,0,0],  # s3 is terminal
         [0,0,0,0,0,0],  # s4 is terminal
         [0,0,0,0,0,0]]
    r = [0, 3, 2, 9] # and denominator 14 to scale back from fractions
    d = 14
    b = calculate_b(m)
    b0 = b[0]
    bd = [round(i*d) for i in b0]
    assert bd == r

    # test 2
    m = [[1,1,1,1], #s0
         [0,0,0,0], #s1
         [0,0,0,0], #s2
         [1,2,3,4]] #s3
    r = [8, 9] # and denominator 17 to scale back from fractions
    d = 17 # denominator to scale back to integers after normalization
    b = calculate_b(m)
    b0 = b[0] # to test transitions only from s0
    bd = [round(i*d) for i in b0]
    assert bd == r
