from Fonctions import *

resultat = bissection("x**2 - 4", 0, 5, 0.001)
print("Racine approchée par la bissection est :", resultat)


# Exemple : chercher la racine de f(x) = x^2 - 4
resultat = balayage("x**2 - 4", 0, 5, 0.1)
print("Racine approchée par balayage est :", resultat)

# Pour Lagrange (Méthode de la sécante)
print("Racine approchée de Lagrange est :", Lagrange("x**2 - 4", 0, 5, 0.001))  # Cherche racine ≈ 2


# Pour Newton
print("Racine approchée de Newton est :", Newton("x**2 - 4", "2*x", 3, 0.0001))  # Cherche racine ≈ 2



"""--------------------pivo de Gauss-----------------"""
A = np.array([[2, -1, 1],
              [3, 3, 9],
              [3, 3, 5]])

b = np.array([8, 0, -6])

solution = gauss_pivot(A, b)
print("Solution :", solution)