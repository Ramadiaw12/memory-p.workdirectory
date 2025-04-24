from Bio import SeqIO  # Importation de Biopython pour lire les fichiers FASTA

motif = "ATG"  # Le motif à rechercher

# Parcours de séquence du fichier FASTA
for record in SeqIO.parse("data/example.fasta", "fasta"):
    seq = str(record.seq)
    positions = []  # Initialisation d'une liste vide pour stocker les positions du motif
    for i in range(len(seq)):
        # On regarde à chaque position de la séquence si les 4 letters qui suivent correspondent au motif
        if seq[i:i+len(motif)] == motif:
            positions.append(i)  # Si oui, on ajoute la position à la liste
        # else:
        #     print("Aucune motif trouvé")
    print(f"{record.id} : {positions}")  # On affiche l'identifiant de la séquence et la liste des positions
