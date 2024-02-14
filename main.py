from methods import *
import numpy as np

def readMatrix():
  rows = []
  with open('./matrix.txt', 'r') as f:
    lines = f.readlines()
    for line in lines:
      numbers = np.array([int(n) for n in line.split()])
      rows.append(numbers)
    f.close()
  return np.mat(rows, float)

def inputGauss(mat: np.matrix):
  try:
    print(gaussElimination(mat))
  except (ZeroDivisionError):
    if mat.item((-1, -1)) != 0 and mat.item((-1, -2)) == 0:
      print('Sistema sem solução')
    else:
      print('Infinitas soluções')

def inputJacobi(mat: np.matrix):
  # read data that you might need, like x(0) and delta
  # put jacobi function as a method in methods.py
  # call jacobi method here after reading inputs
  pass

def inputSeidel(mat: np.matrix):
  # read data that you might need, like x(0) and delta
  # put seidel function as a method in methods.py
  # call seidel method here after reading inputs
  pass

inputs = {
  "1": inputGauss,
  "2": inputJacobi,
  "3": inputSeidel
}

if __name__ == '__main__':
  print("\n*** DIGITE A MATRIZ NO ARQUIVO matrix.txt ***\n")
  print("Escolha o método a ser utilizado")
  print(" 1 - Eliminação de Gauss")
  print(" 2 - Gauss-Jacobi")
  print(" 3 - Gauss-Seidel\n")
  inputs[input("Método: ")](readMatrix())
