import numpy as np

np.set_printoptions(formatter={'float_kind':'{:.4f}'.format})

def gaussElimination(mat: np.matrix):
  rows = mat.shape[0]
  cols = mat.shape[1]
  xs = np.zeros(rows)
  indexXs = [i for i in range(rows)] # if columns are swapped

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
    leftSide = mat.item((i, cols-1))
    for j in range(cols-2, i, -1):
      leftSide -= mat.item(i, j)*xs[indexXs[j]]
    xs[indexXs[i]] = leftSide/mat.item(i, i)

  s = ''
  for i in range(len(xs)):
    s += f'x{i+1}: {xs[i]:.4f}\n'
  return s

# Possível
# mat = np.mat([[1, 2, 10], [2, -1, 5]], float)

# Impossível
# mat = np.mat([[1, 1, 10], [1, 1, 5]], float)

# Infinitas soluções
# mat = np.mat([[1, 1, 3], [2, 2, 6]], float)

# try:
#   print(gaussElimination(mat))
# except (ZeroDivisionError):
#   if mat.item((-1, -1)) != 0 and mat.item((-1, -2)) == 0:
#     print('Sistema sem solução')
#   else:
#     print('Infinitas soluções')