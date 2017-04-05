# absorbing-markov-chainsasdfasdfsadfdsf
Pure Python 2.7 implementation of solving Absorbing Markov Chains (no dependencies)

## Motivation
Matrix operations in pure Python are nothing complex but boring.

I guess you're looking for implementation to run in Python 2.7 sandbox. In this case this set of functions might save you several hours of coding and bugfixing.

This set of functions can do:
* Identity matrix generation
* Matrix substruction
* Matrix multiplications
* Matrix inverse
* Check if the matrix is zero-matrix
* Decompose Markov matrix on Q and R components

These are basically what you need to calculate absorbing probabilities, e.g. `B = (I - Q)^-1 * R` 

## Background
