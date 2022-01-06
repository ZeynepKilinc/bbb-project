import pandas as pd
import glob
import os
from scipy.stats import shapiro
from matplotlib import pyplot
from numpy import mean



aas = ["T","I","W","M","-","Y","R","L","N","C","A","D","E","K","H","G","V","Q","P","F","S"]
def fasta2dict(filename):
    #take fasta as dict
    seqDict = {}
    filein = open(filename, "r")
    for line in filein:
        if line.startswith(">"):
            header = line[1:].strip()
            seqDict[header] = ""
        else:
            seqDict[header] += line.strip()
    return seqDict

def seq2conserv(seqDict):    
    #calculate conservations in a dict
    maindict={}
    for pos in range(len(list(seqDict.values())[0])):
        maindict[pos]=[0]*len(aas)
        for seq in seqDict.keys():
            aa=seqDict[seq][pos]
            if aa in aas:
                i=aas.index(aa)
                maindict[pos][i]+=1
        for sum in range(len(maindict[pos])):
            maindict[pos][sum]=maindict[pos][sum]/len(seqDict.keys())
    return maindict

def conserv2file(filename,motif,protein):
    seqDict=fasta2dict(filename) 
    maindict=seq2conserv(seqDict)
    df = pd.DataFrame(columns = aas)
    for pos in maindict.keys():
        tempdict={}
        for aap in range(len(aas)):
            tempdict[aas[aap]]=maindict[pos][aap]
        df = df.append(tempdict, ignore_index=True)
    df.to_csv("../data/"+motif+"/conservation/"+motif+"_"+protein+"_conservation.tsv", sep="\t",index=False)
    # for pos in maindict.keys():
    #     for aap in range(len(aas)):
    #         aaper.write(str(pos)+ "\t" + str(aas[aap])+"\t"+ str(maindict[pos][aap])+"\n")
    # refseq=seqDict["A0AVI4"]
    # return maindict,refseq
            

def writeall(motif):
    flist=glob.glob("../data/"+motif+"/alignment/*.fasta")
    print("Number of proteins with ortholog files for",motif,"is",len(flist),"unique proteins",len(list(set(flist))))
    os.system("mkdir ../data/"+motif+"/conservation")
    for i in range(0,len(flist),1):
        print(flist[i])
        proteinFile = flist[i]
        protein=proteinFile.split("/")[4][:-10]
        conserv2file(proteinFile,motif,protein)


def getconsv(data):
    #get the highest consv for each position 
    stat, p = shapiro(data)
    print('Statistics=%.3f, p=%.3f' % (stat, p))
    # interpret
    alpha = 0.05
    if p > alpha:
        print('Sample looks Gaussian (fail to reject H0)')
    else:
        print('Sample does not look Gaussian (reject H0)')  
            
def readconservationfile(motif,protein,filename,f):
    df = pd.read_csv("../data/"+motif+"/conservation/"+motif+"_"+protein+"_conservation.tsv", sep='\t')
    consvlist={}
    for i in range(len(df)):
        m=max(df.loc[i,])
        idx=list(df.loc[i,]).index(m)
        if idx!=4: #!= "-" not gap
            consvlist[i]=[m,idx,"aa"]
        else:
            consvlist[i]=[m,idx,"g"]
        
    maxswgaps = [value[0] for key, value in consvlist.items() ]
    maxs = [value[0] for key, value in consvlist.items() if value[2]!="g" ]
    
    #print(consvlist)
    df2 = pd.read_csv("../data/"+motif+"/"+motif+"_proteins.tsv", sep='\t')
    rows=df2[df2['UniproId'].str.contains(protein)]
    seqDict=fasta2dict(filename)
    
    msaseq = [value for key, value in seqDict.items() if key == protein ]
    print(msaseq)
    msaseq=msaseq[0]

    for i in range(len(rows)):
        print(rows["motifstart"].values[i])
        motifstart,motifend,gapcount,aacount,msapos,motifscores=rows["motifstart"].values[i],rows["motifend"].values[i],0,0,[],[]
        for a in range(len(msaseq)):
            if msaseq[a] =="-":
                gapcount+=1
            else:
                aacount+=1
                if motifstart<=aacount<=motifend:
                    msapos.append(gapcount+aacount-1)
                    motifscores.append([consvlist[gapcount+aacount-1][0],msaseq[a]])
        print(msapos,gapcount,aacount)
        msamean,msameanwogaps,motifmean,motifscore=mean(maxswgaps),mean(maxs),0,[]
        for t in range(len(motifscores)):
            if motif.split("-")[t]=="X":
                motifscore.append("X")
            elif motifscores[t][1]==motif.split("-")[t]:
                motifmean+=motifscores[t][0]
                motifscore.append(motifscores[t][0])
            else:
                motifscore.append(-motifscores[t][0])
        motifmean=motifmean/(len(motifscore)-motifscore.count("X"))
        motifper= len([item for item in maxswgaps if item > motifmean ])/len(maxswgaps)*100
        motifperwogaps= len([item for item in maxs if item > motifmean ])/len(maxs)*100
        motifperwogaps2= len([item for item in maxs if item > motifmean ])/len(maxswgaps)*100
        f.write(protein+"\t"+str(msamean)+"\t" +str(msameanwogaps)+
            "\t"+str(motifmean)+"\t"+str(motifscore)+"\t"+str(motifper)+"\t"+str(motifperwogaps)+"\t"+str(motifperwogaps2)+"\n")
            
              
    #stat, p = shapiro(consvlist)
    #pyplot.bar(list(consvlist.keys()), maxs, color ='maroon')
    #pyplot.show()
    # if p>0.05:
    #     print('Sample looks Gaussian (fail to reject H0)')
    # else:
    #     print(p)
      

    #read consv file as df
    #read protein file as df
    #get the lines of protein
    #get the positions of motif



    

# Shapiro-Wilk Test
def runthefolder(motif,f):
    f.write("UniproId"+ "\t"+"MSA_mean"+"\t" +"MSA_mean_wo_gaps"+
        "\t"+"motifMean"+"\t"+"motifScore"+"\t"+"motifPercent%"+"\t"+"motifPercent_wogaps%"+"\t"+"motifPercent_wogaps2%"+"\n")
       
    flist=glob.glob(motif+"_msa/*.fasta")
    for i in range(0,len(flist),1):
        fname=flist[i]
        protein=fname.split("\\")[1].split("_")[0]
        print(fname)
        filename=fname.replace("\\","/")
        readconservationfile(motif,protein,filename,f)
# protein="A1KXE4"
# filename=motif+"_msa/"+protein+"_orthologs_msa.fasta"
# f = open(motif+"_proteins/"+motif.replace("-","")+"_results.tsv", "w")
# f.write("UniproId"+ "\t"+"MSA_mean"+"\t" +"MSA_mean_wo_gaps"+
#         "\t"+"motifMean"+"\t"+"motifScore"+"\t"+"motifPercent%"+"\t"+"motifPercent_wogaps%"+"\t"+"motifPercent_wogaps2%"+"\n")
# readconservationfile(motif,protein,filename,f)
# filein = open("NPXYproteins.tsv", "r")
# lines=[]
# for line in filein:
#     lines.append(line.strip())
# for i in range(1,len(lines)):
#     line=lines[i]
#     if line.split("\t")[4]== "A0AVI4":
#         motifstart=line.split("\t")[0]
#         motifend=line.split("\t")[1]

# gapcount=0
# aacount=0
# for i in range(len(refseq)):
#     if refseq[i] =="-":
#         gapcount+=1
#     else:
#         aacount+=1
#     if aacount==motifstart:
#      print(refseq[i])

    

