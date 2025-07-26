from Bio import SeqIO
import os
import sys

# Recupération du chemin FASTA passé à l'argument au moment de l'éxécution
file_path = sys.argv[1]

file_path = sys.argv[1]
if len(sys.argv) < 2:
    print("Veuillez vérifier le nombre d'argument ecrite. Exemple: python filename.py data/monfichier.fasta")
    sys.exit(1)
else:
    
    if not os.path.isfile(file_path):
        print(f"ERREUR: le fichier {file_path} n'existe pas ou est invalide")
        sys.exit(1)

# Paramètres de filtrage. les séquences entre 200 et 1000
min_length = 50  # minimum 100 nucléotides
max_length = 100  # ou None pour pas de maximum

total_sequences = 0
sequences_conserves = 0

# Lecture des sequences
for record in SeqIO.parse(file_path, "fasta"):

    print(f"La séquence {record.id} a pour longueur {len(record.seq)}")
    # print(f"l'ID de la sequence est {record.id}")
    total_sequences += 1 # Le nombre de séquence que contient le fichier fasta
    seq_length = len(record.seq) # La longueur de chaque séquence
    
    # Filtrage des séquences
    if seq_length >= min_length:  
        if max_length is None or seq_length <= max_length: # Conserve toutes séquences supérieur à min_length et les seq inférieur ou égal à max_length
            sequences_conserves += 1 # Compte les sequences concervées 

# Affichage des résultats
print(f"Nombre total de séquences : {total_sequences}")
print(f"Nombre de séquences conservées après filtrage : {sequences_conserves}")
