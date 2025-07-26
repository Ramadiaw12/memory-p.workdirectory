


import sys
# def rama(nom="Mon nom", age=0):    
#     print(nom, age)
# nom="Mon nom"
# age=0
# if len(sys.argv) > 1:
#     nom=sys.argv[1]
#     age=sys.argv[2]
# rama(nom, age)
#creation de fonction
def filtre_genome(filtre):
    print(f'genomes filtré avec une taille {filtre}')
#Programme principal
filtre=50
if len(sys.argv) > 1:
    if int(sys.argv[1]):
        filtre=int(sys.argv[1])
filtre_genome(filtre)