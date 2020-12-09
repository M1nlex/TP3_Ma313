import numpy as np
import time
import matplotlib.pyplot as plt



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


def MIRelaxation(a, b, x0, epsilon, nitermax, w = 1):
    d = np.diag(np.diag(a))
    e = -np.tril(a - d)
    f = -np.triu(a - d)
    #w = 1
    m = (1/w)*d - e
    n = ((1/w)-1)*d + f
    return migenerale(m, n, b, x0, epsilon, nitermax)


def erreur(a, b, x0, epsilon, nitermax, methode, type_erreur):

    if methode == 1:
        x, i, e = mijacobi(a, b, x0, epsilon, nitermax)
    if methode == 2:
        x, i, e = MIGaussSeidel(a, b, x0, epsilon, nitermax)
    if methode == 3:
        x, i, e = MIRelaxation(a, b, x0, epsilon, nitermax, w=1)
    if type_erreur == 1:
        return np.mean(np.abs(np.dot(a, x) - b))
    if type_erreur == 2:
        return np.linalg.norm(np.dot(a, x) - b)


def erreur_graph(n_max=100, nb_matrice=10, epsilon=10**-6, nitermax=100, methode=2, type_erreur=1):

    e_moyen = []
    x = []
    for n in range(2, n_max):
        print(n)
        liste_e = []
        for j in range(0, nb_matrice):
            # Création matrices aléatoires
            a = diagonale_dominante(n)
            b = np.transpose(np.random.randn(1, n))
            x0 = np.ones((n, 1))

            e = erreur(a, b, x0, epsilon, nitermax, methode, type_erreur)
            liste_e.append(e)
        x.append(n)
        e_moyen.append([np.mean(liste_e)])

    plt.plot(x, e_moyen)
    plt.xlabel("Taille")
    plt.ylabel("Erreur")
    plt.show()



def test1():
    print("oui")


def w_optimal(n_max = 10, epsilon = 10**(-7) , nitermax = 100 , nb_matrix_max = 10 , precision_i = 10 , type = 0):
    for n in range(2,n_max):
        L_time_moyen = []
        L_iterr_moyen = []
        L_vconverg_moyen = []
        L_i = []
        for i in np.linspace(1 , 2 , precision_i):
            L_time = []
            L_iterr = []
            L_vconverg = []
            for j in range(0,nb_matrix_max):
                #Création matrices aléatoires
                A = diagonale_dominante(n)
                B = np.transpose(np.random.randn(1,n))
                x0 = np.ones((n,1))

                #Temps de calcul
                t0 = time.perf_counter()
                x,iterr,eps = MIRelaxation(A, B, x0, epsilon, nitermax , w = i)
                t1 = time.perf_counter()
                t = t1 - t0

                L_time.append(t)
                L_iterr.append(iterr)
                if iterr!=0:
                    L_vconverg.append(1/iterr)
            #Ajout de la valeur de w
            L_i.append(i)

            #Calcul et ajout du temps moyen sur 'nb_matrix_max' matrices
            t_moyen = sum(L_time)/nb_matrix_max
            L_time_moyen.append(t_moyen)

            #Calcul et ajout nbre itération moyen
            iterr_moyen = sum(L_iterr)/nb_matrix_max
            L_iterr_moyen.append(iterr_moyen)

            #Calcul vitesse convergence moyenne
            vconverg_moyen = sum(L_vconverg)/nb_matrix_max
            L_vconverg_moyen.append(vconverg_moyen)

        if type == 0:
            #Calcul du w minimum par le temps
            t1_mini = min(L_time_moyen)
            indice1 = L_time_moyen.index(t1_mini)
            w_opti1 = L_i[indice1]
            print("Le w_opti par le temps pour la taille "+str(n)+" est : "+str(w_opti1))

            #Courbes par le temps
            plt.plot(L_i, L_time_moyen, ".:", label = "Temps pour la taille "+str(n))

        if type == 1:
            #Calcul du w minimum par le temps
            t2_mini = min(L_iterr_moyen)
            indice2 = L_iterr_moyen.index(t2_mini)
            w_opti2 = L_i[indice2]
            print("Le w_opti par le nbr d'itérations pour la taille "+str(n)+" est : "+str(w_opti2))

            #Courbes le nbr d'itérations
            plt.plot(L_i, L_iterr_moyen, ".:", label = "Nbr itér pour la taille "+str(n))

        if type == 2:
            plt.plot(L_i, L_vconverg_moyen, ".:", label = "Vitesse de convergence pour la taille "+str(n))
    if type == 0:
        plt.ylabel("Temps de calcul (s)")
        plt.title("Temps nécessaire pour la résolution d'un système en fonction du coefficient w\n", fontsize=12)
    if type == 1 :
        plt.ylabel("Nbr d'itérations")
        plt.title("Nombre d'itérations pour la résolution d'un système en fonction du coefficient w\n", fontsize=12)
    if type == 2:
        plt.ylabel("Vitesse de convergence")
        plt.title("Vitesse de convergence pour la résolution d'un système en fonction du coefficient w\n", fontsize=12)
    plt.xlabel("valeur de w")
    plt.legend(loc = "upper left")

    plt.show()

def rayon_spectral(A):
    L1, L2 = np.linalg.eig(A)
    R = max( [abs(number) for number in L1] )
    return R

def w_parfait(A):
    #Calcul de M^-1 * N pour le rayon spectral
    d = np.diag(np.diag(A))
    e = -np.tril(A-d)
    f = -np.triu(A-d)
    m = d
    n = e + f
    m1 = np.linalg.inv(m)
    B = np.dot(m1, n)
    #Rayon spectral et w
    r = rayon_spectral(B)
    w = 1 + (( r/( 1+np.sqrt( 1-(r**2) ) ) )**2)
    return w

#Création des matrices : A TESTER

def Matrice_A1 (n=100):
    A = np.zeros((n,n))
    for i in range (n):
        for j in range (n):
            if i != j :
                A[i,j] = 1/(12 + (3*1 - 5*j)**2)
            elif i == j :
                A[i,j] = 3
    print (A)
    return A

def Vecteur_b (n=100):
    b = np.zeros((n,1))
    for i in range (n):
        b[i,0] = np.cos(i/8)
    print (b)
    return b

def Matrice_A2 (n=100):
    A = np.zeros((n,n))
    for i in range (n):
        for j in range (n):
            A[i,j] = 1/(1 + 3*abs(i-j))
    print (A)
    return A

def diagonale_dominante(n,x=1):
    A = np.zeros((n,n))
    for i in range (0,n):
        S = 0
        for j in range (0,n):
            A[i,j] = x*np.random.random()
            S += abs(A[i,j])
        A[i,i] = S + x
    #print (A)
    return A


#Test (Alexandre)
matrice_test_A = np.array([[1,3,2,4],[5,3,0,4],[10,7,2,3],[4,7,8,2]])
matrice_test_X_resu = np.transpose( np.array([[7,8,15,6]]) )
matrice_test_X = np.transpose( np.array([[1,1,1,1]]) )
matrice_test_B = np.transpose( np.array([[85,83,174,216]]) )

#Autre test (Romaric)
matrice_test_A1 = [[2, 1], [5, 7]]
matrice_test_B1 = np.array([[11, 13]]).T
matrice_test_X1 = np.array([[1, 1]]).T
matrice_test_X1_resu = np.array([[7.111, -3.222]]).T

#Encore un autre test (Marianne)
A_test = [[2,1,0],[1,2,1],[0,1,2]]
X_test_resu = np.transpose( np.array([[1,2,3]]) )
X_test = np.transpose( np.array([[1,1,1]]) )
B_test = np.transpose( np.array([[4,8,8]]) )
