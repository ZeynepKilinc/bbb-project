

import requests
import os



def getortho(motif):
    print("Ortholog sequences will be retrieved from OMA database for each protein that has the motif in its cytosolic region and has a transmembrane domain")
    print("----------getortho.py started-----------")
    filein = open("../data/"+motif+"/"+motif+"_proteins.tsv", "r")
    lines=[]
    os.system("mkdir ../data/"+motif+"/orthologs")
    for line in filein:
        lines.append(line.strip())


    for i in range(1,len(lines)):
        line=lines[i]
        unikod=line.split("\t")[4]
        url= "https://omabrowser.org/oma/omagroup/"+unikod+"/fasta/"
        r = requests.get(url)
        #s = open(motif+"_orthologs/"+unikod+".fasta", "w")
        #s.write("> "+unikod+"\n"+line.split("\t")[9])
        #print(line.split("\t")[9])
        #s.close()
        content = r.content
        content=content.decode("utf-8")
        sequences=content.split("\n")
        f = open("../data/"+motif+"/orthologs/"+unikod+"_orthologs.fasta", "w")
        sequ=""
        for seq in sequences:
            if len(seq)>0:
                if seq[0]==">":
                    if sequ!="":
                        f.write(sequ+"\n")
                    f.write(seq+"\n")
                    sequ=""
                else:
                    sequ+=seq.strip()
        f.write(sequ+"\n")
        f.write("> "+unikod+"\n"+line.split("\t")[9])
        f.close()
        if i==int(len(lines)/4):
            print("----%25 completed----")  
        if i==int(len(lines)/2):
            print("----%50 completed----")
        if i==int(len(lines)/4*3):
            print("----%75 completed----") 
        if i==int(len(lines)):
            print("----%100 completed----") 
    print("----------getortho.py ended-------------\n\n")





#url= "https://omabrowser.org/api/protein/Q7Z408/orthologs/"

#use requests



