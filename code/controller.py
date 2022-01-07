import os
from getuniprot import getuniprot
from getortho import getortho
from makealignment import makealignment
from conservation import writeall,runthefolder

#This is the code to run all other python scripts and terminal commands

"""First lets run the getuniprot.py which will first 
    get the proteins that has the searched motif from Expasy in the;
    sp|Q14738|2A5D_HUMAN    516     519     USERPAT1        .       .       .      NPqY
    format.

    Then it will search the uniprot id's. After retrieving the uniprot id 
    it will check if the found motif position is in a cytoplasmic region and
    if the protein has a transmembrane region. If both conditions are provided 
    it will print the;
    uniprotID, ensemlID, orthodID, omaID, genecardsID, motif start position, 
    motif end position, cytopolasmic region start position,cytopolasmic region end position 
    To a tsv file named as {motif}proteins.tsv
"""
if __name__ == "__main__":
    print("######### Pipeline Started ################")
    repeat=int(input("Enter the number of motifs you want to search for:"))
    for rep in range(repeat):
        motif=input("Enter the "+ str(rep+1)+". motif you are searching for :")
        print("#########################")
        getuniprot(motif)
        getortho(motif)
        makealignment(motif)
        writeall(motif)
        f = open("../data/"+motif+"/"+motif+"_results.tsv", "w")
        runthefolder(motif,f)
        print("######### Your "+str(rep+1)+". "+motif+ "'s folders are ready ################")


#clustalo --profile1 orthologs/A0AVI4.fasta -i orthologs/A0AVI4_orthologs.fasta -o A0AVI4_orthologs_out.fasta --force --auto

