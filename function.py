import numpy as np


def migenerale(m, n, b, x0, epsilon, nitermax):
    i = 0
    e = 10
    while i <= nitermax and e > epsilon:
        i += 1
        u = n * x0 + b

        x1 = np.linalg.solve(m, u)

        e = np.linalg.norm(x1-x0)
        x0 = x1

    return x0, i, e


def test1():
    print("oui")



def Matrice_A (n):
    A = np.zeros(n,n)
    for i in range (n):
        for j in range (n):
            if i != j :
                A[i,j] = 1/(12 + (3*1 - 5*j)**2)
            elif i == j :
                A[i,j] = 3
    print (A)
    return A
