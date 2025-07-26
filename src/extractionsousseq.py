from Bio import SeqIO
import os 
import sys
import gzip

# Recupération du chemin FASTq passé à l'argument au moment de l'éxécution
file_path = sys.argv[1]

if len(sys.argv) < 2:
    print("Veuillez vérifier le nombre d'argument ecrite. Exemple: python filename.py data/monfichier.fasta")
    sys.exit(1)

# Vérifie si le fichier existe et est un fichier
if not os.path.isfile(file_path):
    print(f"Erreur : le fichier '{file_path}' n'existe pas ou est invalide.")
    sys.exit(1)

file_path = sys.argv[1]
 
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

start = 10
end = 30

    
sous_sequence = record.seq[start-1:end]  # -1  les index Python commencent à 0
print(f"La sequence {sequence_choisie.id} a pour sous-séquence ({start}-{end}) : {sous_sequence}")
