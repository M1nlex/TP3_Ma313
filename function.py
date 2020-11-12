import numpy as np


def migenerale(m, n, b, x0, epsilon, nitermax):
    i = 0
    e = 10
    while i <= nitermax and e > epsilon:
        i += 1
        u = n * x0 + b

        x1 = np.linalg.solve(m, u)

        e = abs(x1-x0)
        x0 = x1

        print(i)

    return x0, i, e


def test1():
    print("oui")
