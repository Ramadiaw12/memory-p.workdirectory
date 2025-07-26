from Bio import SeqIO
import os
import sys
import gzip

file_path = sys.argv[1]

if len(sys.argv) < 2:
    print("Veuillez vérifier le nombre d'argument ecrite. Exemple: python filename.py data/monfichier.fasta")
    sys.exit(1)

file_path = sys.argv[1]

# Vérifier si le fichier donné par le user est un fasta ou fastq
if file_path.endswith(".fasta.gz"):
    file_format = "fasta"
    handle = gzip.open(file_path, "rt")  # lecture texte
elif file_path.endswith(".fasta"):
    file_format = "fasta"
    handle = open(file_path, "r")  # lecture
elif file_path.endswith(".fastq.gz"):
    file_format = "fastq"
    handle = gzip.open(file_path, "rt")  # lecture texte
elif file_path.endswith(".fastq"):
    file_format = "fastq"
    handle = open(file_path, "r")  # lecture
else:
    print("Format non reconnu. Utilisez un fichier .fasta ou .fastq")
    sys.exit(1) 

# Vérifie si le fichier existe et est un fichier
if not os.path.isfile(file_path):
    print(f"Erreur : le fichier '{file_path}' est introuvable.")
    sys.exit(1)

# Parcours et lit le fichier . Charge toutes les séquences dans une liste
sequences = list(SeqIO.parse(file_path, file_format))
# Affiche la liste des IDs disponibles
print("Séquences disponibles :")
for idx, record in enumerate(sequences):
    print(f"{idx + 1}. ID = {record.id}")

# Demande à l'utilisateur de choisir une séquence
choix = input("Entrez le numéro ou l'ID de la séquence à analyser : ")

sequence_choisie = None

# Si l'utilisateur entre un numéro
if choix.isdigit():
    index = int(choix) - 1
    if 0 <= index < len(sequences):
        sequence_choisie = sequences[index]
    else:
        print("Numéro invalide.")
        sys.exit(1)
else:
    # Si l'utilisateur entre un ID
    for record in sequences:
        if record.id == choix:
            sequence_choisie = record
            break
    if sequence_choisie is None:
        print("ID non trouvé.")
        sys.exit(1)

# Déclaration et Initialisation des bases
base_counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
# Parcourir notre fichier fasta et lit la séquence

for base in sequence_choisie.seq:  # Parcours la séquence (record.seq issue de l'objet Seqreccord)
    if base in base_counts:  # Check si la base est dans base_counts
        base_counts[base] += 1  # Incrémente de 1 le nombre les bases présentent dans base_counts
# L'affichage des rapports de base dans la séquence
print('Le nombre de base A est :', base_counts['A'])
print('Le nombre de base T est :', base_counts['T'])
print('Le nombre de base C est :', base_counts['C'])
print('Le nombre de base G est :', base_counts['G'])

