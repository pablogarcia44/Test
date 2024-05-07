import openmc

library = openmc.data.DataLibrary()


import os



def lister_chemins_fichiers(chemin):
    # Vérifier si le chemin est un dossier existant
    if not os.path.isdir(chemin):
        print("Le chemin spécifié n'est pas un dossier valide.")
        return []

    # Initialiser une liste pour stocker les chemins d'accès aux fichiers
    chemins_fichiers = []

    # Parcourir tous les fichiers et dossiers dans le dossier donné
    for nom in os.listdir(chemin):
        # Construire le chemin complet du fichier
        chemin_fichier = os.path.join(chemin, nom)
        # Vérifier si le chemin correspond à un fichier
        if os.path.isfile(chemin_fichier):
            chemins_fichiers.append(chemin_fichier)

    return chemins_fichiers



# Obtenir le répertoire du fichier Python en cours d'exécution
repertoire_courant = os.path.dirname(__file__)

# Nom du dossier que vous souhaitez atteindre (remplacez 'nom_dossier' par le nom réel du dossier)
nom_dossier = 'jeff-3.3-hdf5'

# Chemin complet du dossier
chemin_dossier = os.path.join(repertoire_courant, nom_dossier)



# Obtenir la liste des chemins d'accès aux fichiers dans le dossier spécifié
chemins_fichiers_dossier = lister_chemins_fichiers(chemin_dossier)

# Afficher les chemins d'accès aux fichiers
print("Chemins d'accès aux fichiers dans le dossier", chemin_dossier, ":")
for chemin_fichier in chemins_fichiers_dossier:
    library.register_file(chemin_fichier)
    print(chemin_fichier)

    
    


library.export_to_xml()





