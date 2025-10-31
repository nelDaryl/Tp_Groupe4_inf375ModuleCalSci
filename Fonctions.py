from sympy import symbols, sympify, lambdify
import numpy as np

"""____________________#--fonction_bissection--#____________________"""
##--------#Cette fonction va trouver une racine de f(x) dans [a, b] avec une précision de e#-----###

def bissection(fonc, a, b, e):
    # Préparation SymPy
    x = symbols('x')
    f = sympify(fonc)              # transforme l'expression textuelle en fonction sympy
    f = lambdify(x, f, 'math')      # transforme en fonction numérique utilisable

    # Ton code original commence ici
    x = a
    y = b
    s = x
    # si f(x) = 0, on vas prendre s = y
    if f(x) == 0:
        s = x

    while (y - x) >= e and (f(x) * f(y)) < 0:
        s = (x + y) / 2
        if f(x) * f(s) < 0:      # Correction mineure nécessaire ici (f(y) → f(s))
            y = s
        else:
            x = s
    return s


"""____________________#--fonction balayage#--#____________________"""
"""#----##cette methode vas trouver une racine approcher de f(x) dans [a, b] avec un pas de e#-----3#"""
def balayage(fonc, a, b, e):
    """
    fonc : fonction en chaîne de caractères ex: "x**2 - 4"
    a, b  : bornes de balayage
    e     : pas de balayage
    """
    x = symbols('x')
    f = sympify(fonc)               # Conversion en expression SymPy
    f_num = lambdify(x, f, 'math')   # Conversion en fonction numérique Python

    # Début du balayage
    X = a

    # Si f(a) = 0 ⇒ a est racine
    if f_num(X) == 0:
        return X

    # Balayage de l'intervalle
    while (X + e) <= b and (f_num(X) * f_num(X + e)) > 0:
        X += e

    # Si aucun changement de signe, racine introuvable
    if (X + e) > b:
        raise ValueError("Pas de changement de signe !!!")

    # Approximation de la racine
    s = X + e/2
    return s


"""____________________#--Methode_de_Lagrange--#____________________"""
##---------------cette mathode est utiliser pour trouver la  racine d'une fonction f(x) = 0.---------------------##

# ---------- Méthode de Lagrange (Sécante) ----------
def Lagrange(f_str, a, b, e):
    x = symbols('x')
    f = sympify(f_str)
    f = lambdify(x, f, 'math')

    x0 = a
    x1 = b
    s = x0
    if f(x1) == 0:
        s = b
    else:
        while (f(x0) * f(x1) < 0) and abs(x0 - x1) >= e:
            s = (x0 * f(x1) - x1 * f(x0)) / (f(x1) - f(x0))
            if f(s) * f(x1) < 0:
                x0 = s
            else:
                x1 = s
    return s


"""____________________#--Methode_de_Newton--#____________________"""
def Newton(f_str, df_str, x0, e, max_iterat=100):
    x = symbols('x')
    f = sympify(f_str)
    df = sympify(df_str)

    f = lambdify(x, f, 'math')
    df = lambdify(x, df, 'math')

    if f(x0) == 0:
        s = x0
    else:
        if df(x0) != 0:
            s = x0 - (f(x0) / df(x0))
            iteration = 0

            while abs(s - x0) >= e and iteration < max_iterat:
                x0 = s
                if df(x0) != 0:
                    s = x0 - (f(x0) / df(x0))
                    iteration += 1
                else:
                    print("La dérivée est nulle en x0")
                    break
            return s
        else:
            print("La dérivée est nulle au point initial.")
            return None


"""____________________#--Methode_gauss_pivot--#____________________"""
def gauss_pivot(A, b):
    # Convertit A et b en flottants pour éviter les divisions entre entiers
    A = A.astype(float)
    b = b.astype(float)

    # Nombre d'équations (ou de lignes de A)
    n = len(b)

    for i in range(n):
        # Étape 1 : Pivot partiel
        # Cherche la ligne avec le plus grand coefficient dans la colonne i
        max_row = np.argmax(abs(A[i:, i])) + i

        # Permute la ligne i avec la ligne max_row
        A[[i, max_row]] = A[[max_row, i]]
        b[i], b[max_row] = b[max_row], b[i]

        # Étape 2 : Élimination pour créer des zéros sous la diagonale
        for j in range(i + 1, n):
            # Coefficient pour annuler A[j][i]
            coef = A[j][i] / A[i][i]

            # Applique l'élimination sur la ligne j
            A[j] = A[j] - coef * A[i]
            b[j] = b[j] - coef * b[i]


    x = np.zeros(n)  # vecteur des solutions initialisé à zéro

    # On résout les équations en partant de la dernière (n-1) jusqu'à la première (0)
    for i in range(n - 1, -1, -1):
        # Calcule la valeur de x[i]
        x[i] = (b[i] - np.dot(A[i][i+1:], x[i+1:])) / A[i][i]

    return x



