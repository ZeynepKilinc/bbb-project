
from Bio.Align.Applications import ClustalOmegaCommandline 
import glob
import os

def makealignment(motif):
    print("Alignments will be build. The code will retrieve each fasta file in [motif]_orthologs folder and prints the output msa to the [motif]_msa folder")
    print("--------makealignment.py starts----------")
    flist=glob.glob("../data/"+motif+"/orthologs/*.fasta")
    print("Number of proteins with ortholog files for",motif,"is",len(flist),"unique proteins",len(list(set(flist))))
    os.system("mkdir ../data/"+motif+"/alignment")
    for i in range(0,len(flist),1):
        fname=flist[i]
        print(fname)
        fname=fname.split("/")[4][:-16]
        
        os.system("clustalo -i "+"../data/"+motif+"/orthologs/"+fname+"_orthologs.fasta -o "+"../data/"+motif+"/alignment/"+fname+"_msa.fasta --force --auto ")
        #os.system("mafft -i orthologs/"+profile+"_orthologs.fasta > msa/"+profile+"_orthologs_out.fasta ")
        #inform the user about the process
        if i==int(len(flist)/4):
                print("----%25 completed----")  
        if i==int(len(flist)/2):
            print("----%50 completed----")
        if i==int(len(flist)/4*3):
            print("----%75 completed----") 
        if i==int(len(flist)):
            print("----%100 completed----") 
    print("--------makealignment.py enden----------\n\n")



"""
in_file = "orthologs/A0AVI4_orthologs.fasta"
out_file = "A0AVI4_orthologs_out.fasta"
clustalomega_cline = ClustalOmegaCommandline(infile = in_file, outfile = out_file, verbose = True, auto = True)
print(clustalomega_cline)
"""