from Bio import SeqIO
import os
import sys
from Bio.SeqUtils import gc_fraction
import gzip


file_path = sys.argv[1]
if len(sys.argv) < 2:
    print("Veuillez vérifier le nombre d'argument ecrite. Exemple: python filename.py data/monfichier.fasta")
    sys.exit(1)

if not os.path.isfile(file_path):
    print(f"ERREUR: le fichier {file_path} n'existe pas ou est invalide")
    sys.exit(1)
   
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
    
# Parcours et lit le fichier . Charge toutes les séquences dans une liste
sequences = list(SeqIO.parse(handle, file_format))

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

seq = str(sequence_choisie.seq).upper()
gc_percent = gc_fraction(seq) * 100
lenght = len(sequence_choisie)
for record in sequences:
    # taille total des lectures
    total_length = len(record.seq)
    
#     # Le teneur en GC
#     seq = record.seq
#     seq_ids.append(record.id)
#     gc_percent = gc_fraction(seq) * 100
#     gc_per_seq.append(gc_percent)


# for seq_id, gc in zip(seq_ids, gc_per_seq):
#     print(f"GC content of {seq_id} est : {gc:.2f}%")


# Affichage des résultats
print("\n--- Résultats ---")
print(f"ID de la séquence : {sequence_choisie.id}")
# Affiche de la taille des lectures
print(f"La longueur de la sequence choisie est {lenght} bases")
# Le pourcentage de GC
print(f"Le contenu en GC est : {gc_percent:.2f}%")
# Et affiche la longueur total
print(f"Le nombre total de lecture est : {len(sequences)} séquences")

