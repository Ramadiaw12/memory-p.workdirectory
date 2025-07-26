from Bio import SeqIO  # Importation de Biopython pour lire les fichiers FASTA
import sys
import os
import gzip

# motif = "ATG"  # Le motif à rechercher
file_path = sys.argv[1]
motif = sys.argv[2]

if len(sys.argv) != 3:
    print("Erreur : usage correct → python filename.py <fichier.fasta> <motif>")
    sys.exit(1)
# Vérifie si le fichier existe et est un fichier
if not os.path.isfile(file_path):
    print(f"Erreur : le fichier '{file_path}' est introuvable.")
    sys.exit(1)

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



# Parcours et lit le fichier . Charge toutes les séquences dans une liste
sequences = list(SeqIO.parse(handle, file_format))

choice = input("Voulez-vous chercher le motif dans tout le fichier ou uniquement dans une séquence spécifique ? (Make your choice to 'all' or 'specificm_sequence') : ")

# Recherche de motif dans tt le fichier
if choice.lower() == 'all':
    
    # lancer la recherche sur tout le fichier fasta/fastq
    results_by_sequence = {}
    # Résultat
    for record in sequences:
        seq_id = record.id
        current_sequence_string = str(record.seq) # Convertir l'objet Seq en chaîne de caractères

        positions_in_current_sequence = []
        motif_len = len(motif)

        # Rechercher le motif dans la séquence actuelle
        for i in range(len(current_sequence_string) - motif_len + 1):
            if current_sequence_string[i:i + motif_len] == motif:
                positions_in_current_sequence.append(i)
        
        # Si le motif a été trouvé dans cette séquence, ajouter ses positions au dictionnaire
        if positions_in_current_sequence:
            results_by_sequence[seq_id] = positions_in_current_sequence
    
    # Affichage des résultats pour toutes les séquences
    if results_by_sequence:
        print(f"\n--- Résultats de la recherche du motif '{motif}' ---")
        for seq_id, positions_list in results_by_sequence.items():
            print(f"Séquence '{seq_id}': positions motif = [ {', '.join(map(str, positions_list))} ]")
    else:
        print(f"Le motif '{motif}' n'a été trouvé dans aucune séquence du fichier.")

    
# rechercher le motif dans la séquence précise
elif choice.lower() == 'specific_sequence':
    # motif = input("Veuillez saisir la séquence spécifique à rechercher : ")

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
    seq = str(sequence_choisie.seq)
    positions = []
    for i in range(len(seq) - len(motif) + 1):
        # On regarde à chaque position de la séquence si les 4 letters qui suivent correspondent au motif
        # On si la position de la séquence à partir de i correspond au motif
        if seq[i:i+len(motif)] == motif:
            positions.append(i)  # Si oui, on ajoute la position à la liste
    # Résultat
    if positions:
        print(f"Le motif '{motif}' de la sequence {sequence_choisie.id} est trouvé aux positions : {positions}")
    else:
        print(f"Aucun motif '{motif}' trouvé dans la séquence '{sequence_choisie.id}'.")
    # print(f"La liste des positions de la séquences {record.id} est : {positions}")  # On affiche l'identifiant de la séquence et la liste des positions

else:
    print("Choix non reconnu, veuillez saisir 'tout' ou 'sequence'.")