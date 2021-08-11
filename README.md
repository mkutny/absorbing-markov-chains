[![Build Status](https://travis-ci.org/mkutny/absorbing-markov-chains.svg?branch=master)](https://travis-ci.org/mkutny/absorbing-markov-chains)

# Absorbing Markov Chains
Pure Python 2.7 implementation of solving Absorbing Markov Chains (no dependencies)

## Motivation
Matrix operations in pure Python are nothing complex but boring. I saw a lot of code snippets in gists and stackexchange questions but I believe that absence of a solid package is a shame.

I guess you're looking for implementation to run in Python 2.7 sandbox. In this case this set of utility functions might save you several hours of coding and bugfixing.

This set of functions can do:
* Identity matrix generation
* Matrix subtraction
* Matrix multiplications
* Matrix inverse
* Matrix reordering (to swap states to bring matrix to canonical form)
* Check if the matrix is zero-matrix
* Decompose Markov matrix on Q and R components

These are basically what you need to calculate some properties like absorbing probabilities.

## Canonical form
Let an absorbing Markov chain with transition matrix `P` have `t` transient states and `r` absorbing states. Then
```
    [ Q R ]
P = [ 0 I ]
```
where `Q` is square `t`-by-`t` matrix, `R` is `t`-by-`r` matrix, `0` is zero-matrix and `I` is identity matrix.

### Absorbing probabilities
If you linearly transform `P` such as `Q` becomes a zero-matrix then what will appear in place of `R` are final absorbing probabilities.
In other words absorbing probabilities `B` could be calculated as `B = (I - Q)^-1 * R`.

## Notes
* These set of utility functions assume that `I` is a zero-matrix
