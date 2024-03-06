import itertools
import numpy as np

matriz = [
    [1, 2, 3, 9],
    [4, 5, 6, 9],
    [7, 8, 9, 9]
]


for p in itertools.permutations(matriz):
    print(p, "\n")


