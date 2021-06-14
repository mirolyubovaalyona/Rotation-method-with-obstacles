
import numpy as np
import sys
import copy
from numpy import array,identity,diagonal
from math import sqrt





a= [[2, 1, 1],
  [1, 2.5, 1],
   [1, 1, 3] ]
a=np.array(a)
b=[3, 5, 3]
n=len(b)

def Jacobi(A):
    n     = A.shape[0]            # matrix size #columns = #lines
    maxit = 100                   # maximum number of iterations
    eps   = 1.0e-15               # accuracy goal
    pi    = np.pi        
    info  = 0                     # return flag
    ev    = np.zeros(n,float)     # initialize eigenvalues
    U     = np.zeros((n,n),float) # initialize eigenvector
    for i in range(0,n): U[i,i] = 1.0

    for t in range(0,maxit):
         s = 0;    # compute sum of off-diagonal elements in A(i,j)
         for i in range(0,n): s = s + np.sum(np.abs(A[i,(i+1):n]))
         if (s < eps): # diagonal form reached
              info = t
              for i in range(0,n):ev[i] = A[i,i]
              break
         else:
              limit = s/(n*(n-1)/2.0)       # average value of off-diagonal elements
              for i in range(0,n-1):       # loop over lines of matrix
                   for j in range(i+1,n):  # loop over columns of matrix
                       if (np.abs(A[i,j]) > limit):      # determine (ij) such that |A(i,j)| larger than average 
                                                         # value of off-diagonal elements
                           denom = A[i,i] - A[j,j]       # denominator of Eq. (3.61)
                           if (np.abs(denom) < eps): phi = pi/4         # Eq. (3.62)
                           else: phi = 0.5*np.arctan(2.0*A[i,j]/denom)  # Eq. (3.61)
                           si = np.sin(phi)
                           co = np.cos(phi)
                           for k in range(i+1,j):
                               store  = A[i,k]
                               A[i,k] = A[i,k]*co + A[k,j]*si  # Eq. (3.56) 
                               A[k,j] = A[k,j]*co - store *si  # Eq. (3.57) 
                           for k in range(j+1,n):
                               store  = A[i,k]
                               A[i,k] = A[i,k]*co + A[j,k]*si  # Eq. (3.56) 
                               A[j,k] = A[j,k]*co - store *si  # Eq. (3.57) 
                           for k in range(0,i):
                               store  = A[k,i]
                               A[k,i] = A[k,i]*co + A[k,j]*si
                               A[k,j] = A[k,j]*co - store *si
                           store = A[i,i]
                           A[i,i] = A[i,i]*co*co + 2.0*A[i,j]*co*si +A[j,j]*si*si  # Eq. (3.58)
                           A[j,j] = A[j,j]*co*co - 2.0*A[i,j]*co*si +store *si*si  # Eq. (3.59)
                           A[i,j] = 0.0                                            # Eq. (3.60)
                           for k in range(0,n):
                                store  = U[k,j]
                                U[k,j] = U[k,j]*co - U[k,i]*si  # Eq. (3.66)
                                U[k,i] = U[k,i]*co + store *si  # Eq. (3.67)
         info = -t # in case no convergence is reached set info to a negative value "-t"
    return ev,U,t

xx=np.linalg.solve(a, b)

ev,U,t=Jacobi(a)
x=[0]*n
print('Теперь известные элементы: ')
for i in range (n):
        x[i]=b[i]/U[i][i]
print(U)
print(xx)