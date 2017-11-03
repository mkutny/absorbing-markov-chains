# PyTest script
# run with: 'py.test -s -v test_amc.py'
import unittest
from fractions import Fraction
import amc

class AmcTests(unittest.TestCase):
    def test_fraction(self):
        m = [[0, 2, 1, 0, 0],
             [0, 0, 0, 3, 4],
             [0, 0, 0, 0, 0],
             [0, 0, 0, 0, 0],
             [0, 0, 0, 0, 0]]

        n = amc.normalize(m, use_fractions=True)
        self.assertTrue( n[0][1] == Fraction(2,3))
        self.assertTrue( n[1][3] == Fraction(3,7))

        b = AmcTests.calculate_b(m, use_fractions=True)
        b0 = b[0]
        self.assertEqual( Fraction(1,3), b0[0])
        self.assertEqual( Fraction(2,7), b0[1])
        self.assertEqual( Fraction(8,21), b0[2])

    def test_sort(self):
        m = [[0,0],[0,0]]
        r = [[0,0],[0,0]]
        self.assertTrue(amc.sort(m) == r)

        m = [[0,0], [1,2]]
        r = [[2,1], [0,0]]
        self.assertTrue(amc.sort(m) == r)

        m = [[1,2,3],
             [0,0,0],
             [3,2,1]]
        r = [[1,3,2],
             [3,1,2],
             [0,0,0]]
        self.assertTrue(amc.sort(m) == r)

        m = [[1,1,1,1], #s0
             [0,0,0,0], #s1
             [0,0,0,0], #s2
             [1,2,3,4]] #s3
        r = [[1,1,1,1], #s0
             [1,4,2,3], #s3
             [0,0,0,0], #s1
             [0,0,0,0]] #s2
        self.assertTrue(amc.sort(m) == r)

    @staticmethod
    def calculate_b(m,use_fractions=False):
        # B = (I-Q)^-1 * R
        m = amc.sort(m)
        n = amc.normalize(m,use_fractions=use_fractions)
        (q, r) = amc.decompose(n)
        i = amc.identity(len(q))
        s = amc.subtract(i, q)
        v = amc.getMatrixInverse(s)
        b = amc.multiply(v, r)
        return b

    def test_b(self):
        # test 1
        m = [[0,1,0,0,0,1],  # s0, the initial state, goes to s1 and s5 with equal probability
             [4,0,0,3,2,0],  # s1 can become s0, s3, or s4, but with different probabilities
             [0,0,0,0,0,0],  # s2 is terminal, and unreachable (never observed in practice)
             [0,0,0,0,0,0],  # s3 is terminal
             [0,0,0,0,0,0],  # s4 is terminal
             [0,0,0,0,0,0]]
        r = [0, 3, 2, 9] # and denominator 14 to scale back from fractions
        d = 14
        b = AmcTests.calculate_b(m)
        b0 = b[0]
        bd = [round(i*d) for i in b0]
        self.assertTrue( bd == r )

        # test 2
        m = [[1,1,1,1], #s0
             [0,0,0,0], #s1
             [0,0,0,0], #s2
             [1,2,3,4]] #s3
        r = [8, 9] # and denominator 17 to scale back from fractions
        d = 17 # denominator to scale back to integers after normalization
        b = AmcTests.calculate_b(m)
        b0 = b[0] # to test transitions only from s0
        bd = [round(i*d) for i in b0]
        self.assertTrue(bd == r)

    @staticmethod
    def lcm(a, b):
        if a > b:
            greater = a
        else:
            greater = b

        while True:
            if greater % a == 0 and greater % b == 0:
                lcm = greater
                break
            greater += 1

        return lcm

    @staticmethod
    def get_lcm_for(l):
        return reduce(lambda x, y: AmcTests.lcm(x, y), l)

    @staticmethod
    def convert_to_lcd(probs):
        ret = []

        least_common_multiple = AmcTests.get_lcm_for([f.denominator for f in probs])
        for f in probs:
            if f.denominator != least_common_multiple:
                ret.append(Fraction(least_common_multiple / f.denominator * f.numerator, least_common_multiple ) )
            else:
                ret.append(Fraction(f.numerator, least_common_multiple ) )
        return ret



    @staticmethod
    def markov_probabilities(m ):
        probs = AmcTests.calculate_b(m, use_fractions=True)[0]
        return AmcTests.convert_to_lcd(probs)


    def test_others(self):
        #a bunch of random tests I used from another project


        m = [[0, 2, 1, 0, 0], [0, 0, 0, 3, 4], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]]
        self.assertEqual([Fraction(7,21), Fraction(6,21), Fraction(8,21)], AmcTests.markov_probabilities(m))

        m = [[0, 1, 0, 0, 0, 1], [4, 0, 0, 3, 2, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0],
             [0, 0, 0, 0, 0, 0]]
        self.assertEqual([0, Fraction(3,14), Fraction(2,14), Fraction(9,14)], AmcTests.markov_probabilities(m))

        m = [
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        ]
        self.assertEqual([Fraction(1,5), Fraction(1,5), Fraction(1,5), Fraction(1,5), Fraction(1,5)], AmcTests.markov_probabilities(m))

        m = [
            [1, 1, 1, 0, 1, 0, 1, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, 1, 1, 1, 0, 1, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, 1, 0, 1, 1, 1, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, 1, 0, 1, 0, 1, 1, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, 1, 0, 1, 0, 1, 0, 1, 1],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        ]
        self.assertEqual([Fraction(1,3), Fraction(1,6), Fraction(1,6), Fraction(1,6), Fraction(1, 6)],
                         AmcTests.markov_probabilities(m))

        self.assertEqual([Fraction(6,100), Fraction(44,100),
                          Fraction(4,100), Fraction(11,100),
                          Fraction(22,100), Fraction(13, 100)],
                         AmcTests.markov_probabilities([
                             [0, 86, 61, 189, 0, 18, 12, 33, 66, 39],
                             [0, 0, 2, 0, 0, 1, 0, 0, 0, 0],
                             [15, 187, 0, 0, 18, 23, 0, 0, 0, 0],
                             [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                         ]
                         ))

        self.assertEqual([Fraction(1,5), Fraction(1,5), Fraction(1,5), Fraction(2, 5)],
                         AmcTests.markov_probabilities([
                             [0, 0, 0, 0, 3, 5, 0, 0, 0, 2],
                             [0, 0, 4, 0, 0, 0, 1, 0, 0, 0],
                             [0, 0, 0, 4, 4, 0, 0, 0, 1, 1],
                             [13, 0, 0, 0, 0, 0, 2, 0, 0, 0],
                             [0, 1, 8, 7, 0, 0, 0, 1, 3, 0],
                             [1, 7, 0, 0, 0, 0, 0, 2, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                         ]
                         )
                         )


