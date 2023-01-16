import pandas as pd
import glob
import os


def fasta2dict(filename):
    #take fasta as dict
    seqDict = {}
    filein = open(filename, "r")
    for line in filein:
        if line.startswith(">"):
            header = line[:].strip()
            seqDict[header] = ""
        else:
            seqDict[header] += line.strip()
    return [seqDict,header]


def gapdeletion(seqdict,postodelete):
    for key,value in seqdict.items():
        for pos in postodelete:
            value=value[:pos] + "." + value[pos+1:]
        seqdict[key]=value.replace(".","")
    return(seqdict)

def removegaps(motif):
    flist=glob.glob("data/"+motif+"/alignment/*.fasta")
    print("Number of proteins with ortholog files for",motif,"is",len(flist),"unique proteins",len(list(set(flist))))
    os.system("mkdir data/"+motif+"/gapsremoved")
    for i in range(0,len(flist),1):
        proteinFile = flist[i]
        seqdict,ref=fasta2dict(proteinFile)
        refseq=seqdict[ref]
        postodelete=([pos for pos, char in enumerate(refseq) if char == "-"])
        seqdict=gapdeletion(seqdict,postodelete)
        f = open("data/"+motif+"/gapsremoved/"+ref[1:]+"_wogaps_msa.fasta", "w")
        for key,value in seqdict.items():
            f.write(key)
            f.write("\n")
            f.write(value)
            f.write("\n")
        f.close()



        
      
        
        

