import numpy as np
from math import inf, fabs

np.set_printoptions(formatter={'float_kind':'{:.4f}'.format})

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

def criterioConvergencia(A):
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
  
  lista = []
  for posicao in resultado[1]:
    lista.append(A[posicao].tolist()[0])

  matriz = np.mat(lista)
  return (True, matriz)

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
  
  resultadoCC = criterioConvergencia(A)

  if resultadoCC[0]:
    A = resultadoCC[1]
    while erroAbsoluto(x0, x1, it) > delta:
      x0 = x1.copy()
      x1 = calculo(A,b,x0)
      for elemento in x1:
        print(f"{elemento:.4f}", end=' ')
      print()
      it += 1
    
    return ""
  else :
    return("O critério de convergência não foi atendido.")
