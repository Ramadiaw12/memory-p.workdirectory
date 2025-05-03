from Bio import SeqIO  # Importation de Biopython pour lire les fichiers FASTA
import sys
import os
# motif = "ATG"  # Le motif à rechercher
if (len(sys.argv) != 3):
    print("Erreur: veuillez entrer deux arguments un fichier fasta et le motif à check")
    sys.exit(1)  # Stop the program if we aren't enough program
else:
    # file_path = 'chemin/vers/votre/fichier.fasta'

    file_path = sys.argv[1]

    # Vérifie si le fichier existe
    if not os.path.exists(file_path):
        print(f"Erreur : le fichier '{file_path}' n'existe pas.")
        sys.exit(1)

    # Recupération des arguments
    fasta_file = sys.argv[1]  # Premier argument après le nom du fichier: nom du fasta file
    motif = sys.argv[2]  # Deuxième ar
    # Parcours de séquence du fichier FASTA
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq)
        positions = []  # Initialisation d'une liste vide pour stocker les positions du motif
        for i in range(len(seq)):
            # On regarde à chaque position de la séquence si les 4 letters qui suivent correspondent au motif
            # On si la position de la séquence à partir de i correspond au motif
            if seq[i:i+len(motif)] == motif:
                positions.append(i)  # Si oui, on ajoute la position à la liste
            # else:
            #     print("Aucune motif trouvé")
        print(f"{record.id} : {positions}")  # On affiche l'identifiant de la séquence et la liste des positions
