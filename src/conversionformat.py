from Bio import SeqIO
import os 
import sys
import gzip

# Récupération du fichier en argument
input_file = sys.argv[1]

# Vérifie que l'utilisateur a bien passé un argument
if len(sys.argv) < 2:
    print("Veuillez vérifier le nombre d'argument ecrite. Exemple: python filename.py data/monfichier.fasta")
    sys.exit(1)

# Vérifie si le fichier existe et est un fichier
if not os.path.isfile(input_file):
    print(f"Erreur : le fichier '{input_file}' n'existe pas ou est invalide.")
    sys.exit(1)

# Verification si le fichier entré est un fichier FASTQ
if input_file.endswith(".fastq.gz"):
    file_format = "fastq"
    handle = gzip.open(input_file, "rt")  # lecture texte
elif input_file.endswith(".fastq"):
    file_format = "fastq"
    handle = open(input_file, "r")  # lecture
else:
    print("Format non reconnu. Utilisez un fichier .fasta ou .fastq")
    sys.exit(1)
      
# Parcours et lit le fichier . Charge toutes les séquences dans une liste
sequences = list(SeqIO.parse(handle, file_format))

# Création du nom du fichier de sortie
output_file = os.path.splitext(input_file)[0] + "_converted.fasta"

# Conversion
count = SeqIO.convert(input_file, "fastq", output_file, "fasta")
print(f"{count} séquences ont été converties et enregistrées dans '{output_file}'.")

