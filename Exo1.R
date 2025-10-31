# Création du tableau de contingence
boisson <- matrix(c(30, 10, 10,
                    20, 25, 5),
                  nrow = 2,
                  byrow = TRUE)

# Nommer les lignes et les colonnes
rownames(boisson) <- c("Homme", "Femme")
colnames(boisson) <- c("Café", "Thé", "Jus")

# Afficher le tableau
print(boisson)

# Réaliser le test du Khi-deux d'indépendance
test_chi2 <- chisq.test(boisson)

# Afficher les résultats du test
print(test_chi2)

# afficher les effectifs attendus :
print(test_chi2$expected)

# afficher la p-valeur seule :
print(test_chi2$p.value)

