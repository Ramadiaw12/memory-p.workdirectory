from Bio import SeqIO

#  Lit la séquence depuis le fichier.fasta
for record in SeqIO.parse("data/example.fasta", "fasta"):
    print(f"ID de la sequence : {record.id}")  # Affiche l'identifiant de la séquence genre l'id
    print(f"la longueur de la sequence est : {len(record.seq)}")  # Donne la longueur de séquence et il a 'record.seq' comme objet