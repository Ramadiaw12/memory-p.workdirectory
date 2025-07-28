import sys
import os
import gzip
from Bio.SeqUtils import gc_fraction
from Bio import SeqIO

# Détection du format et ouverture du fichier en fonction

def detect_format_and_open(file_path):
    if file_path.endswith(".fasta.gz"):
        return "fasta", gzip.open(file_path, "rt")
    elif file_path.endswith(".fasta"):
        return "fasta", open(file_path, "r")
    elif file_path.endswith(".fastq.gz"):
        return "fastq", gzip.open(file_path, "rt")
    elif file_path.endswith(".fastq"):
        return "fastq", open(file_path, "r")
    else:
        raise ValueError("Format non reconnu. Utilisez un fichier .fasta ou .fastq")



def calcul_gc(file_path):
    file_format, handle = detect_format_and_open(file_path)
    # if len(sys.argv) < 2:
    #     print("Veuillez vérifier le nombre d'arguments. Exemple: python filename.py data/monfichier.fasta")
    #     sys.exit(1)

    # file_path = sys.argv[1]

    # if not os.path.exists(file_path):
    #     print(f"ERREUR: le fichier {file_path} n'existe pas ou est invalide")
    #     sys.exit(1)
    
    
    # Initialisation des variables
    total_gc = 0
    total_length = 0
    nb_sequences = 0

    # Parcours le fichier et lit le fichier
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

# Recherche de motif

# Affichage des positions de motif du fichier
def find_motif_all(sequences, motif):
    results = {}
    motif = motif.upper() # Les motifs en lettre majuscule
    motif_len = len(motif)
    for record in sequences:
        seq_id = record.id # L'ID de la sequence
        sequence_str = str(record.seq).upper() # C'est la sequence
        positions = [ # Création d'une liste de position 
            i for i in range(len(sequence_str) - motif_len + 1) # Génération de toutes les motifs possible oû la motif longueur peut commencer dans la sequence
            if sequence_str[i:i + motif_len] == motif # Check dans la sous-chaîne(motif_len) à partir de la position i est égale au motif rechercher
        ]
        if positions:
            results[seq_id] = positions
        # Affichage UNE SEULE FOIS après avoir traité toutes les séquences
    if results:
        print(f"\nMotif '{motif}' trouvé dans {len(results)} séquence(s) :")
        for idx, (seq_id, positions) in enumerate(results.items(), start=1):
            print(f"Les positions de la sequence {seq_id} sont  positions : {positions}")
    else:
        print(f"\n Aucune séquence ne contient le motif '{motif}'.")

    return results
    

# Affichage des motifs de la sequence choisie
def find_motif_in_sequence(sequences, motif, target_id):
    for record in sequences:
        if record.id == target_id:
            seq = str(record.seq)
            # Creation d'une dictionnaire pour afficher la position des motifs
            positions = [
                i for i in range(len(seq) - len(motif) + 1)
                if seq[i:i + len(motif)] == motif
            ]
            return record.id, positions
    raise ValueError("ID de séquence non trouvé.")

def search_motif(file_path, motif, mode='all', specific_id=None):
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"Fichier '{file_path}' introuvable.")

    file_format, handle = detect_format_and_open(file_path)
    sequences = list(SeqIO.parse(handle, file_format))
    handle.close()

    if mode == 'all':
        result = find_motif_all(sequences, motif)
        return result
    # Recherche de motif dans la sequence specifique
    elif mode == 'specific_sequence':
        if specific_id is None:
            raise ValueError("Vous devez fournir l’ID de la séquence pour ce mode.")
        seq_id, positions = find_motif_in_sequence(sequences, motif, specific_id)
        
        if positions:
            print(f"\nMotif '{motif}' trouvé dans la séquence {seq_id} aux positions : {positions}")
        else:
            print(f"\nMotif '{motif}' NON trouvé dans la séquence {seq_id}.")
        
        return {seq_id: positions}
    
# Affichage de la taille et du nombre de sequence
def size_nb_seq (file_path):
    file_format, handle = detect_format_and_open(file_path)
    # Initialisation des variables
    total_len = 0  
    nb_seq = 0

    for record in SeqIO.parse(handle, file_format):
        seq = record.seq
        nb_seq +=1
        total_len += len(seq)
    print(f"Nombre total de séquences : {nb_seq}")
    print(f"Taille totale des séquences : {total_len} bases")

# Choix de l'utilisateur
# La fonction principale
def main():
    # Vérification des fichiers entrés
    if len(sys.argv) < 2:
        print("Erreur : Veuillez vérifier le nombre d'arguments. Exemple: chemin/vers/fichier.fasta(.gz)")
        sys.exit(1)

    file_path = sys.argv[1]
    if not os.path.exists(file_path):
        print(f"ERREUR: le fichier {file_path} n'existe pas ou est invalide")
        sys.exit(1)

    print("\nQue souhaitez-vous faire ?")
    print("1. Calculer la teneur en GC")
    print("2. Rechercher un motif dans toutes les séquences")
    print("3. Rechercher un motif dans une séquence spécifique")
    print("4. Rechercher la taille et le nombre de sequence")


    choix = input("Votre choix (1/2/3/4) : ")

    # A ffichage du GC content
    if choix == "1":
        calcul_gc(file_path)
    # Affichage des positions du motif du fichier
    elif choix == "2":
        motif = input("Entrez le motif à rechercher : ").upper()
        search_motif(file_path, motif, mode='all')
    # Afiichage des positions du motif de la sequence choisit
    elif choix == "3":
        motif = input("Entrez le motif à rechercher : ").upper()
    


        # Lire les séquences
        file_format, handle = detect_format_and_open(file_path)
        sequences = list(SeqIO.parse(handle, file_format))

        # Vérifier qu'on a des séquences
        if not sequences:
            print("Aucune séquence trouvée dans le fichier.")
            sys.exit(1)

        # Afficher les IDs numérotés
        print("\nSéquences disponibles :")
        for idx, record in enumerate(sequences, start=1):
            print(f"{idx}. {record.id}")

        # Choisir via un numéro
        choix_num = input("\nEntrez le numéro de la séquence ciblée : ").strip()
        while not choix_num.isdigit() or int(choix_num) < 1 or int(choix_num) > len(sequences):
            print("Entrée invalide. Veuillez entrer un numéro valide.")
            choix_num = input("Entrez le numéro de la séquence ciblée : ").strip()

        # Obtenir l’ID correspondant au numéro choisi
        specific_id = sequences[int(choix_num) - 1].id

        # Lancer la recherche
        search_motif(file_path, motif, mode='specific_sequence', specific_id=specific_id)
    # Affichage de la taille et du nombre de sequence
    elif choix == "4":
        size_nb_seq(file_path)
    else:
        print("Choix invalide.")

if __name__ == "__main__":
    main()