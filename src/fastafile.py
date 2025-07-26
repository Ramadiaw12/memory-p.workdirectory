import sys
import os
import gzip
from Bio.SeqUtils import gc_fraction
from Bio import SeqIO

def statistic():
    if len(sys.argv) < 2:
        print("Veuillez vérifier le nombre d'arguments. Exemple: python filename.py data/monfichier.fasta")
        sys.exit(1)

    file_path = sys.argv[1]

    if not os.path.isfile(file_path):
        print(f"ERREUR: le fichier {file_path} n'existe pas ou est invalide")
        sys.exit(1)
    
    # Détection du format et ouverture du fichier en fonction
    if sys.argv[1].endswith(".fasta.gz"):
        file_format = "fasta"
        handle = gzip.open(sys.argv[1], "rt")  # lecture texte
    elif sys.argv[1].endswith(".fasta"):
        file_format = "fasta"
        handle = open(sys.argv[1], "r")  # lecture
    else:
        print("Format non reconnu. Utilisez un fichier .fasta ou .fasta.gz")
        sys.exit(1)  
    # Initialisation des variables
    total_gc = 0
    total_length = 0
    nb_sequences = 0

    for record in SeqIO.parse(handle, file_format):
        seq = record.seq
        gc = gc_fraction(seq)
        gc_percent = gc * 100

        # Le pourcentage de GC de chaque sequence
        print(f"Le pourcentage de GC de {record.id} est de : {gc_percent:.2f}%")

        # Le teneur en GC global
        total_gc += gc * len(seq)
        total_length += len(seq)
        nb_sequences += 1

    if nb_sequences == 0:
        print("Aucune séquence trouvée dans le fichier.")
        sys.exit(1)

    global_gc_percent = (total_gc / total_length) * 100

    print(f"Nombre de séquences analysées : {nb_sequences}")
    print(f"Teneur globale en GC : {global_gc_percent:.2f}%")
        
if __name__ == "__main__":
    statistic()