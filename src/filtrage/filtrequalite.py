from Bio import SeqIO
import os
import sys



# Recupération du chemin FASTQ passé à l'argument au moment de l'éxécution
file_path = sys.argv[1]
# Verifie l'existane du fastq
file_path = sys.argv[1]
if len(sys.argv) < 2:
    print("Veuillez vérifier le nombre d'argument ecrite. Exemple: python filename.py data/monfichier.fasta")
    sys.exit(1)
else:
    
    if not os.path.isfile(file_path):
        print(f"ERREUR: le fichier {file_path} n'existe pas ou est invalide")
        sys.exit(1)

min_quality = 20
min_length = 30
for record in SeqIO.parse(file_path, "fastq"):
    mySeq = record.seq
    # Calcule de la qualité moyenne
    avg_quality = sum(record.letter_annotations["phred_quality"]) / len(record)
    
    # Crontrole de séquence si il passe les controle
    if avg_quality >= min_quality and len(record) >= min_length:
        print(f"Séquence {record.id} retenue. La qualité moyenne est : {avg_quality:.1f}")   
    
    else:
        print(f"Séquence {record.id} ignorée (qualité = {avg_quality:.1f}, longueur = {len(record)})")