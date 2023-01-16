import requests
import os
import pandas as pd

def getpdb(motif):
    df,pdblink = pd.read_csv("../data/"+motif+"/"+motif+"_proteins.tsv", sep='\t'),[]
    for i in range(len(df["UniproId"])):
        url= "https://alphafold.ebi.ac.uk/files/AF-"+df.at[i,"UniproId"]+"-F1-model_v4.pdb"
        pdblink.append(url)
        r = requests.get(url)
        print(i)
        #open("P01130.pdb", "wb").write(r.content)
    df.insert(8, "pdbLink", pdblink)
    df.to_csv("../data/"+motif+"/"+motif+"_proteins.tsv", sep="\t",index=False)
   
getpdb("N-P-X-Y")