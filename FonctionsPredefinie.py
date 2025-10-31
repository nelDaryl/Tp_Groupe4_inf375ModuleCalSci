from sympy import symbols, sympify, lambdify

"""----------------------------BIssection------------------------"""
def bissection(f_str, a, b, e):
    x = symbols('x')
    f = lambdify(x, sympify(f_str))

    x0 = a
    x1 = b
    s = x0

    if f(x0) == 0:
        s = x0
    else:
        while (x1 - x0) >= e and f(x0) * f(x1) < 0:
            s = (x0 + x1) / 2
            if f(x0) * f(s) < 0:
                x1 = s
            else:
                x0 = s
    return s



"""----------------------------Lagrange------------------------"""
def Lagrange(f_str, a, b, e):
    x = symbols('x')
    f = lambdify(x, sympify(f_str))

    x0 = a
    x1 = b
    s = x0

    if f(x1) == 0:
        s = b
    else:
        while f(x0)*f(x1) < 0 and abs(x0 - x1) >= e:
            s = (x0*f(x1) - x1*f(x0)) / (f(x1) - f(x0))
            if f(s)*f(x1) < 0:
                x0 = s
            else:
                x1 = s
    return s



"""----------------------------Newton------------------------"""
def Newton(f_str, df_str, x0, e, max_iterat=100):
    x = symbols('x')
    f = lambdify(x, sympify(f_str))
    df = lambdify(x, sympify(df_str))

    if f(x0) == 0:
        return x0
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
                    print("⚠️ La dérivée est nulle en x0.")
                    break
            return s
        else:
            print("⚠️ La dérivée est nulle en x0.")
            return None


"""----------------------------Pivo de Gauss------------------------"""
import numpy as np

def gauss_pivot(A, b):
    A = np.array(A, float)
    b = np.array(b, float)
    n = len(b)

    # Élimination
    for k in range(n - 1):
        pivot = np.argmax(abs(A[k:n, k])) + k
        if A[pivot, k] == 0:
            raise ValueError("Pivot nul — pas de solution unique")
        if pivot != k:
            A[[k, pivot]] = A[[pivot, k]]
            b[[k, pivot]] = b[[pivot, k]]
        for i in range(k + 1, n):
            factor = A[i, k] / A[k, k]
            A[i, k:] -= factor * A[k, k:]
            b[i] -= factor * b[k]

    # Substitution arrière
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x[i] = (b[i] - np.dot(A[i, i+1:], x[i+1:])) / A[i, i]
    return x




"""----------------------------Balayage------------------------"""
def balayage(f_str, a, b, e):
    """
    Méthode du balayage pour trouver une racine de f(x)=0 sur [a, b]
    ---------------------------------------------------------------
    f_str : str -> expression de la fonction (ex: "x**2 - 4")
    a, b : bornes de l'intervalle
    e : pas ou seuil de précision
    """
    x = symbols('x')
    f = lambdify(x, sympify(f_str))  # Conversion en fonction Python

    X = a
    i = 0
    s = None  # Racine approximée

    # Si la fonction est nulle à gauche
    if f(X) == 0:
        s = X
    else:
        # Tant qu'on reste dans l'intervalle et que le signe ne change pas
        while (X + e) <= b and (f(X) * f(X + e)) > 0:
            X = X + e
            i += 1

        # On prend la racine approximative au milieu du dernier intervalle
        s = X + e / 2

    return s



"""--------------------------test sur les modules--------------------------------"""

result = bissection("x**2 - 4", 0, 3, 1e-6)
print("Racine ≈", result)

result = Lagrange("x**3 - 2*x - 5", 2, 3, 1e-6)
print("Racine ≈", result)


result = Newton("x**3 - 2*x - 5", "3*x**2 - 2", 2, 1e-6)
print("Racine ≈", result)


A = [[2, 1, -1],
     [-3, -1, 2],
     [-2, 1, 2]]
b = [8, -11, -3]

result = gauss_pivot(A, b)
print("Solution :", result)

result = balayage("x**3 - 2*x - 5", 0, 5, 0.1)
print("Racine approximative ≈", result)
