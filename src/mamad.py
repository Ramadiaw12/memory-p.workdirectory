import sys
import os
from Bio.SeqUtils import gc_fraction



def teneur_GC(file, seq="*"):
    """
    python3 def_arg.py seq

    seq "*"     Permet d'appliquer le pourcentage de GC dans une sequence specifique du fichier ou dans outes les sequences (e.g. ID)
    """
    
    #calcul du GC total
    DNA=str(file)
  
    GC = (DNA.count('G') + DNA.count('C')) / len(DNA)

    Percent_GC = GC * 100
    return round(Percent_GC, 2)


if __name__=="__main__":

    print(teneur_GC.__doc__)
    if sys.argv[2]=="ID":
        file=sys.argv[1]

