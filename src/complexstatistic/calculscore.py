from Bio import SeqIO
import os
import sys


file_path = sys.argv[1]

if not os.path.exists(file_path):
    print(f"ERREUR: le fichier {file_path} n'est pas un fichier valide")
sys.exit(1)