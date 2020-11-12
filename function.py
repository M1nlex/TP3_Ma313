import numpy as np


def migenerale(m, n, b, x0, epsilon, nitermax):
    i = 0
    e = 10
    while i <= nitermax and e > epsilon:
        i += 1
        u = n @ x0 + b
        #print(u)

        x1 = np.linalg.solve(m, u)

        e = np.linalg.norm(x1-x0)
        x0 = x1

    return x0, i, e


def mijacobi(a, b, x0, epsilon, nitermax):

    m = np.diag(np.diag(a))
    n = m - a
    #print(a, m, n)

    print(migenerale(m, n, b, x0, epsilon, nitermax))






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

matrice_test_A = [[1,3,2,4],[5,3,0,4],[10,7,2,3],[4,7,8,2]]
matrice_test_X = [[7,8,15,6]]
matrice_test_B = [[85,83,174,216]]
