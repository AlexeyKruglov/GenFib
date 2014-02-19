GenFib
======

Generalized Fibonacci sequence fast calulator

Calculate N first elements of a generalized Fiblonacci sequence defined by
  a_0 = 1
  a_i = sum_{k=1}^i w_k a_{i-k}
in Theta(N log^2 N) time and Theta(N) memory. N is assumed to be a power of 2 in the code.
