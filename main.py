from methods import *
import numpy as np

def readMatrix():
  rows = []
  with open('./matrix.txt', 'r') as f:
    lines = f.readlines()
    for line in lines:
      numbers = np.array([float(n) for n in line.split()])
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
  # Adjusting data read to fit function format
  A = []
  b = []
  for i in range(mat.shape[0]):
    lista = []
    for j in range(mat.shape[1]):
      if j == mat.shape[1] - 1:
        b.append(mat[i,j])
      else:
        lista.append(mat[i,j])
    A.append(lista)
  A = np.mat(A)
  b = np.array(b)

  x0 = np.array([])
  # Reading data needed
  print("Insira os valores do vetor x0:")
  x1 = []
  for i in range(b.shape[0]):
    x1.append(float(input()))
  x1 = np.array(x1)

  print()
  delta = float(input("Insira a precisão desejada: "))
  
  print(jacobi(A, b, x0, x1, delta))

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

def main():
  print("\n*** DIGITE A MATRIZ NO ARQUIVO matrix.txt ***\n")
  print("Escolha o método a ser utilizado")
  print(" 1 - Eliminação de Gauss")
  print(" 2 - Gauss-Jacobi")
  print(" 3 - Gauss-Seidel\n")
  inputs[input("Método: ")](readMatrix())

if __name__ == '__main__':
  main()
