from typing import List
from math import inf, fabs
import numpy as np
import itertools

np.set_printoptions(formatter={'float_kind':'{:.4f}'.format})

# *** Método da Eliminação Gaussiana ***

def gaussElimination(mat: np.matrix):
   rows = mat.shape[0]
   cols = mat.shape[1]
   xs = np.zeros(cols-1)
   indexXs = [i for i in range(cols-1)] # if columns are swapped

   for i in range(rows-1):
      if mat.item(i,i) == 0:
         # switch rows
         for j in range(i+1, rows):
            if mat.item(j,i) != 0:
               mat[[i, j]] = mat [[j, i]]
               break
         if mat.item(i,i) == 0:
            # switch cols if couldnt switch rows
            for j in range(i, cols-1):
               if mat.item(i,j) != 0:
                  mat[:, i], mat[:, j] = mat[:, j], mat[:, i].copy()
                  indexXs[i] = j
                  indexXs[j] = i
                  break
      if mat.item(i,i) != 0:
         # upper triangular matrix
         for j in range(i+1, rows):
            mj = mat.item(j,i)/mat.item(i,i)
            mat[j] = mat[j]-mj*mat[i]
            mat.itemset((j, i), 0)

   for i in range(rows-1, -1, -1):
      rightSide = mat.item((i, len(xs)))
      idx = (cols-2)-(rows-1-i)
      for j in range(cols-2, idx, -1):
         rightSide -= mat.item(i, j)*xs[indexXs[j]]
      xs[indexXs[idx]] = rightSide/mat.item(i, idx)

   s = ''
   for i in range(len(xs)):
      s += f'x{i+1}: {xs[i]:.4f}\n'
   return s

# *** Método Gauss-Jacobi ***

def verificaFormato(A, b):
   return A.shape[0] == A.shape[1] and A.shape[0] == len(b)

def verificacaoRecursiva(possiveisPosicoes, indice, vetor):
   if indice < len(possiveisPosicoes) - 1:
      for posicao in possiveisPosicoes[indice]:
         if posicao not in vetor:
            vetor.append(posicao)
            resultado = verificacaoRecursiva(possiveisPosicoes, indice + 1, vetor.copy())
            if resultado:
               return resultado
            vetor.remove(posicao)
      return (False, None)
   else:
      for posicao in possiveisPosicoes[indice]:
         if posicao not in vetor:
            vetor.append(posicao)
            return (True, vetor)
      return (False, None)

def criterioConvergencia(A, b : np.array):
   possiveisPosicoes = []

   for i in range(A.shape[0]):                              # Para cada linha da matriz
      posicoes = []                                          # Verifica as posições que ela pode ocupar
      for j in range(A.shape[1]):                            # Para cada elemento que pode estar na diagonal principal
         valor = 0                                            # Zera o valor do cálculo
         for k in range(A.shape[1]):                          # Para cada elemento que não vai estar na diagonal principal
            if j != k:
               try:
                  valor += fabs(A[i,k]) / fabs(A[i,j]) # Faz o cálculo
               except:
                  valor += inf    
         if valor < 1:
            posicoes.append(j)
      if len(posicoes) == 0 or (len(posicoes) == 1 and posicoes in possiveisPosicoes):
         return(False, None)
      possiveisPosicoes.append(posicoes)

   resultado = verificacaoRecursiva(possiveisPosicoes, 0, [])

   if not resultado[0]:
      return (False, None)
  
   # lista = [None]
   # for posicao in resultado[1]:
   #    lista.append(A[posicao].tolist()[0])

   lista = [None for i in range(A.shape[0])]
   b = b.tolist()
   bNovo = b.copy()
   for i in range(len(resultado[1])):
      lista[resultado[1][i]] = A[i].tolist()[0]
      bNovo[resultado[1][i]] = b[i]

   matriz = np.mat(lista)
   bNovo = np.array(bNovo)
   return (True, matriz, bNovo)

def calculo(A, b, x0):
   x1 = []
   for i in range(A.shape[0]):
      valor = 0
      for j in range(A.shape[1]):
         if i != j:
            valor += x0[j] * A[i,j]
      try:
         x1.append((b[i] - valor)/(A[i,i]))
      except:
         print("Erro: divisão por 0")
         exit()
   return x1

def erroAbsoluto(x0, x1, it):
   if it == 0:
      return inf
   maiorNum = 0
   maiorDen = 0
   for i in range(len(x0)):
      dif = fabs(x1[i] - x0[i])
      if dif > maiorNum:
         maiorNum = dif
   for i in range(len(x1)):
      if fabs(x1[i]) > maiorDen:
         maiorDen = fabs(x1[i])
   try:
      erro = maiorNum / maiorDen
   except:
      erro = inf
   return erro

def jacobi(A : np.matrix, b : np.array, x0 : np.array, x1 : np.array, delta : float):
   if(not verificaFormato(A, b)):
      return("Erro de formato") 
   it = 0
   
   resultadoCC = criterioConvergencia(A, b)

   if resultadoCC[0]:
      A = resultadoCC[1]
      b = resultadoCC[2]
      while erroAbsoluto(x0, x1, it) > delta:
         x0 = x1.copy()
         x1 = calculo(A,b,x0)
         for elemento in x1:
            print(f"{elemento:.4f}", end=' ')
         print()
         it += 1
         
      return ""
   else:
      return("O critério de convergência não foi atendido.")

# *** Método Gauss-Seidel *** 

def gaussSeidel(mat: np.matrix, delta: float, x_k: np.array) -> str:
   numberOfRows, numberOfColumns = mat.shape
   if not gaussSeidelConvergenceCriterion(mat, numberOfRows):
      mat = gaussSeidelSwappRows(mat)
      if mat is None:
         return "Não é possível solucionar o sistema linear através do método Gauss Seidel"
   
   x_k = gaussSeidelInitialAproximation(mat, x_k)

   x_kPlus1 = np.array([0 for z in range(numberOfRows)], float)
   while True:
      x_kPlus1 = np.array([0 for z in range(numberOfRows)], float)
      for x in range(numberOfRows):
         x_kPlus1[x] = mat.item(x, numberOfColumns-1)
         for var in range(numberOfRows): 
            if var != x:
               x_kPlus1[x] -= mat.item(x, var) * (x_kPlus1[var] if var < x else x_k[var])
         x_kPlus1[x] /= mat.item(x, x)
      if (max(np.abs((x_k - x_kPlus1) / (max(np.abs(x_kPlus1))))) <= delta): break
      x_k = x_kPlus1.copy()
   
   #return x_kPlus1
   
   stringifiedResult = ''
   for i in range(len(x_kPlus1)):
      stringifiedResult += f'x{i+1}: {x_kPlus1[i]:.4f}\n'
   return stringifiedResult

def gaussSeidelInitialAproximation(mat: np.matrix, x_k: np.array) -> np.array:
   if x_k.size == 0:
      x_k = np.zeros(mat.shape[0], float)
      for i in range(mat.shape[0]):
         x_k[i] = (mat.item(i, -1) / mat.item(i, i))
   else: 
      return x_k

def gaussSeidelConvergenceCriterion(mat: np.matrix, size: int) -> bool:
   beta = np.zeros(size, float)
   for i in range(size):
      for j in range(size): 
         if i != j:
            beta[i] += (abs(mat.item(i, j)) * (beta[j] if j < i else 1))
      if abs(mat.item(i, i)) == 0:
         beta[i] = 2
      else:
         beta[i] /= abs(mat.item(i, i))
   return max(beta) < 1

def gaussSeidelSwappRows(mat: np.matrix) -> List:
   numberOfRows, numberOfColumns = mat.shape
   matAuxiliar = mat.copy()

   for p in itertools.permutations(mat):
      newMat = []
      for line in p:
         l = []
         for i in range(line.shape[1]):
            l.append(line.item(0,i))
         newMat.append(l)
      newMat = np.mat(newMat)
      if gaussSeidelConvergenceCriterion(newMat, numberOfRows):
         return np.matrix(newMat)
   
   return None