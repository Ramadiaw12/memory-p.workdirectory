from Bio import SeqIO
import sys
import os
import gzip
from Bio.SeqUtils import gc_fraction


file_path = sys.argv[1]
# Check le fichier est ce que c'est un fichier valide
if not os.path.exists(file_path):
    print(f"Erreur : le fichier '{file_path}' n'est pas un fichier valid")
    sys.exit(1)
else:
    with gzip.open(file_path, "rt", encoding="utf-8") as handle:  # Ouvre le fichier .gz . handle: c'est le fichier ouvert ouvert que nous allons lire
        # Parcours la séquence et la lit une par une. enumarate compte la séquence en commençant par 1 avec 'i' le numéro de la séquence
        for i, record in enumerate(SeqIO.parse(handle, "fastq"), 1):
            # print(record.seq)
            print(f"Séquence {i} : longueur = {len(record.seq)}")  # Affiche le numéro de la séquence et sa la longueur
            # Affiche le tenneur de GC dans la séquence
            print(f"GC content of {record.seq} : {gc_fraction(record.seq):.2f}")

# Initialisé le compteur à 0
total = 0
total += len(record.seq)
print(f"la longueur total  de la séquence est : {total}")
