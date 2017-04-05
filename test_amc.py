# PyTest script
# run with: 'py.test -s -v test_amc.py'

import pytest
import amc

def test_sort():
    m = [[1,2,3],
         [0,0,0],
         [3,2,1]]
    r = [[1,3,2],
         [3,1,2],
         [0,0,0]]
    assert amc.sort(m) == r

    m = [[0,0], [1,2]]
    r = [[2,1], [0,0]]
    assert amc.sort(m) == r

    m = [[0,0],[0,0]]
    r = [[0,0],[0,0]]
    assert amc.sort(m) == r
    
    m = [[1,1,1,1],
         [0,0,0,0],
         [0,0,0,0],
         [1,2,3,4]]
    r = [[1,1,1,1],
         [1,4,2,3],
         [0,0,0,0],
         [0,0,0,0]]
    assert amc.sort(m) == r
