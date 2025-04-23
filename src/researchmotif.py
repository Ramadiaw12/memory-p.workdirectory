from Bio import SeqIO

motif = "TATA"
for record in SeqIO.parse("example.fasta", "fasta"):
    print(record.id)