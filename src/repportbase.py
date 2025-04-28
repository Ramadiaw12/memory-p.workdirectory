from Bio import SeqIO

# Déclaration et Initialisation des bases
base_counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
# Parcourir notre fichier fasta et lit la séquence
for record in SeqIO.parse("data/example.fasta", "fasta"):
    for base in record.seq:  # Parcours la séquence (record.seq issue de l'objet Seqreccord)
        if base in base_counts:  # Check si la base est dans base_counts
            base_counts[base] += 1  # Incrémente de 1 le nombre les bases présentent dans base_counts
# L'affichage des rapports de base dans la séquence
print('Le nombre de base A est :', base_counts['A'])
print('Le nombre de base T est :', base_counts['T'])
print('Le nombre de base C est :', base_counts['C'])
print('Le nombre de base G est :', base_counts['G'])

