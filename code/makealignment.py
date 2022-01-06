
from Bio.Align.Applications import ClustalOmegaCommandline 
import glob
import os

def makealignment(motif):
    print("Alignments will be build. The code will retrieve each fasta file in [motif]_orthologs folder and prints the output msa to the [motif]_msa folder")
    print("--------makealignment.py starts----------")
    flist=glob.glob(motif+"_orthologs/*.fasta")
    print("Number of proteins with ortholog files for",motif,"is",len(flist),"unique proteins",len(list(set(flist))))
    os.system("mkdir "+motif+"_msa")
    for i in range(0,len(flist),1):
        fname=flist[i]
        fname=fname.split("/")[1][:-16]
        print(fname)
        os.system("clustalo -i "+motif+"_orthologs/"+fname+"_orthologs.fasta -o "+motif+"_msa/"+fname+"_orthologs_msa.fasta --force --auto ")
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