from Bio import SeqIO

# Valid sequences
seq_nucleotide = set("ATCG")

for record in SeqIO.parse("data/example.fasta", "fasta-blast"):
    # Regarder dans la séquences
    sequence = str(record.seq).upper()

# Check if the sequences
    if set(sequence).issubset(seq_nucleotide):
        print(f"{record.id} : Séquence VALIDE")
    else:
        print(f"{record.id} : Séquence NON VALIDE")
