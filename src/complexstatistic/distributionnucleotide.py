from Bio import SeqIO
import os
import sys
from collections import defaultdict, Counter

file_path = sys.argv[1]

if len(sys.argv[1]) < 2:
    print("Veuillez vérifier le nombre d'argument ecrite. Exemple: python filename.py data/monfichier.fasta")
    sys.exit(1)

else:
    if not os.path.exists(file_path):
        print(f"ERREUR: le fichier {file_path} n'est pas un fchier valide")
        sys.exit(1)

if file_path.endswith(".fasta"):
    file_format = "fasta"
elif file_path.endswith(".fastq"):
    file_format = "fastq"
else:
    print("Format non reconnu. Utilisez un fichier .fasta ou .fastq")
    sys.exit(1)

# Parcours et lit le fichier . Charge toutes les séquences dans une liste
sequences = list(SeqIO.parse(file_path, file_format))

# Initialisation : dictionnaire {position: liste des bases}
position_nucleotides = defaultdict(list)

# Remplir le dictionnaire position_nucleotides
for record in sequences:
    seq = str(record.seq).upper()
    for i, base in enumerate(seq):
        position_nucleotides[i].append(base)

# Analyse de la distribution
print("Distribution des nucléotides par position :\n")
for position in sorted(position_nucleotides.keys()):
    base_counts = Counter(position_nucleotides[position])
    total = sum(base_counts.values())
    # print(f"Séquence {record.id}")
    print(f"Position {position + 1} : ", end="")
    for base in ['A', 'T', 'G', 'C']:
        count = base_counts.get(base, 0)
        percent = (count / total) * 100 if total > 0 else 0
        print(f"{base}={count} ({percent:.1f}%)", end="")
