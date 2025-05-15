from Bio import SeqIO
import os
import sys
import gzip

file_path = sys.argv[1]

if not os.path.exists(file_path):
    print(f"Erreur : le fichier '{file_path}' n'est pas un fichier valid ")
    sys.exit(1)
else:
    # Valid sequences
    for record in SeqIO.parse(gzip.open(file_path, 'rt', encoding='utf-8'), "fastq"):
        print("C'est bien un fichier fastq valid!")

# for record in SeqIO.parse(file_path2, 'fasta'):
#     print("C'est bien un fichier fasta valid !!!")
