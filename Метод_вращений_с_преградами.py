
import numpy as np
import sys
import copy




a= [[2, 1, 1],
  [1, 2.5, 1],
   [1, 1, 3] ]
a=np.array(a)
b=[3, 5, 3]
n=len(b)

def f(n, a, b):
    x=[0]*n
    for k in range(2):
        for l in range(n):
            for i in range(n):
                for j in range(n):
                    d=(((a[i][i]-a[j][j])**2)+4*(a[i][j]**2))**(1/2)
                    c=(1/2*(1+(abs(a[i][i]-a[j][j]))/d))**(1/2)
                    sgn=np.sign(a[i][j]*(a[i][i]-a[j][j]))
                    s=sgn*((1/2*(1-(abs(a[i][i]-a[j][j]))/d))**(1/2))
                    a[i][i]=(c**2)*a[i][i]+2*c*s*a[i][j]+(s**2)*a[j][j]
                    a[j][j]=(s**2)*a[i][i]-2*c*s*a[i][j]+(c**2)*a[j][j]
                    a[i][j]=a[j][i]= ((s**2)-(c**2))*a[i][j]+s*c*(a[i][i]-a[j][j])
                    
                    if (k!=i) and (k!=j):
                        a[k][i], a[i][k]=c*a[k][i]+s*a[k][j], c*a[k][i]+s*a[k][j]
                        a[k][j], a[j][k]=-s*a[k][i]+c*a[k][j], -s*a[k][i]+c*a[k][j]
                    if (k!=i) and (k!=j) and (l!=i) and (l!=j):
                        a[k][l]=a[k][l]

                    
                    
        
    print(a.round(8))
    
    for i in range(n):
       x[i]=b[i]/a[i][i]
        
    return x

xx=np.linalg.solve(a, b)

x=f(n, a, b)
# Displaying solution
print('Теперь известные элементы: ')
print(x)
print(xx)