import numpy as np


def migenerale(m, n, b, x0, epsilon, nitermax):
    i = 0
    e = 10
    while i <= nitermax and e > epsilon:
        i += 1
        u = np.dot(n, x0) + b
        x1 = np.linalg.solve(m, u)

        e = np.linalg.norm(x1-x0)
        x0 = x1

    return x0, i, e


def mijacobi(a, b, x0, epsilon, nitermax):
    d = np.diag(np.diag(a))
    e = -np.tril(a-d)
    f = -np.triu(a-d)
    m = d
    n = e + f
    return migenerale(m, n, b, x0, epsilon, nitermax)


def MIGaussSeidel(a, b, x0, epsilon, nitermax):
    d = np.diag(np.diag(a))
    e = -np.tril(a-d)
    f = -np.triu(a-d)
    m = d - e
    n = f
    return migenerale(m, n, b, x0, epsilon, nitermax)


def MIRelaxation(a, b, x0, epsilon, nitermax):
    d = np.diag(np.diag(a))
    e = np.tril(a)-d
    f = np.triu(a)-d
    w = 1
    m = (1/w)*d - e
    n = ((1/w)-1)*d + f
    return(migenerale(m, n, b, x0, epsilon, nitermax))

def test1():
    print("oui")



#Création des matrices : A TESTER

def Matrice_A1 (n):
    A = np.zeros(n,n)
    for i in range (n):
        for j in range (n):
            if i != j :
                A[i,j] = 1/(12 + (3*1 - 5*j)**2)
            elif i == j :
                A[i,j] = 3
    print (A)
    return A

def Vecteur_b (n):
    b = np.zeros(n,1)
    for i in range (n):
        b[i,1] = cos(i/8)
    print (b)
    return b

def Matrice_A2 (n):
    A = np.zeros(n,n)
    for i in range (n):
        for j in range (n):
            A[i,j] = 1/(1 + 3*abs(i-j))
    print (A)
    return A


#Test :

matrice_test_A = np.array([[1,3,2,4],[5,3,0,4],[10,7,2,3],[4,7,8,2]])
matrice_test_X = np.transpose( np.array([[7,8,15,6]]) )
matrice_test_B = np.transpose( np.array([[85,83,174,216]]) )

#Autre test
matrice_test_A1 = [[2, 1], [5, 7]]
matrice_test_B1 = np.array([[11, 13]]).T
matrice_test_X1 = np.array([[1, 1]]).T

A_test = [[2,1,0],[1,2,1],[0,1,2]]
X_test = np.transpose([[1,2,3]])
B_test = np.transpose([[4,8,8]])
