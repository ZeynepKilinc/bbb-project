import requests
import os


def getuniprot(motif):
    print("Pipeline started with motif "+ motif)
    print("The code getuniprot.py will be execuded, there will be one output file '[motif]_proteins.tsv'")
    print("---------getuniprot.py started---------")
    os.system("mkdir ../data/"+motif)
    f = open("../data/"+motif+"/"+motif+"_proteins.tsv", "w")
    f.write("motifstart"+ "\t" +"motifend"+ "\t" +"siteStart"+ "\t" +"siteEnd"+"\t"+
                "UniproId"+ "\t" +"EnsemblId"+"\t" +"OrthodbId"+"\t"+"GeneCardsId"+"\t"+"OmaId"+"\t"+"Sequence"+"\n")

    url= "https://prosite.expasy.org/cgi-bin/prosite/PSScan.cgi?sig="+motif+"&lineage=Homo%20sapiens&db=sp&output=txt"
    #use requests
    r = requests.get(url)
    content = r.content
    content=content.decode("utf-8") 
    content1=content.splitlines() 
    #orthodb
    #alternatif olarak ensembl
    for l in range(0,len(content1),1):
        line1=content1[l]
        hastransmem=0
        iscytop=0
        finalrange,ensemblid,orthodb,genecards,omaId=[],"","","",""
        if len(line1)>0:
            #sp|Q14738|2A5D_HUMAN    516     519     USERPAT1        .       .       .      NPqY
            print(line1)
            head=line1.split("\t")[0].split("|")[1]
            url= "https://www.uniprot.org/uniprot/"+head+".txt"
            r = requests.get(url)
            motifstart=line1.split("\t")[1]
            motifend=line1.split("\t")[2]
            content = r.content
            content=content.decode("utf-8") 
            content=content.splitlines() 
            sequencestartl=100**2
            for t in range(len(content)):
                line=content[t]
                line=line.split("   ")
                if line[0] =="FT":
                    if line[1]=="TRANSMEM":
                        hastransmem=1
                    if line[1]=="TOPO_DOM":
                        if "Cytoplasmic" in content[t+1].split("   ")[6]:
                            rng=line[3].replace(" ","").split("..")
                            if int(rng[0])<int(motifstart) and int(motifend)<int(rng[1]):
                                iscytop=1
                                finalrange=rng
                if line[0] =="DR":
                    if line[1].split(";")[0] =="Ensembl":
                        ensemblid=line[1].split(";")[2]
                if line[0] =="DR":
                    if line[1].split(";")[0] =="GeneCards":
                        genecards=line[1].split(";")[1]
                if line[0] =="DR":
                    if line[1].split(";")[0] =="OrthoDB":
                        orthodb=line[1].split(";")[1]
                if line[0] =="DR":
                    if line[1].split(";")[0] =="OMA":
                        omaId=line[1].split(";")[1]
                if line[0] =="SQ":
                    sequence=""
                    sequencestartl=t
                if t>sequencestartl and line[0]!="//":
                    sequence+=content[t].strip()
                    sequence=sequence.replace(" ","")
            if hastransmem==1 and iscytop==1:
                #print(len(finalrange),motifstart,motifend,head,finalrange[0]+ "\t" +finalrange[1],ensemblid)
                f.write(motifstart+ "\t" +motifend+ "\t" +finalrange[0]+ "\t" +finalrange[1]+"\t"+
                            head+"\t"+ensemblid+"\t"+orthodb+"\t"+genecards+"\t"+omaId+"\t"+sequence+"\n")
            if l==int(len(content1)/4):
                print("----%25 completed----")  
            if l==int(len(content1)/2):
                print("----%50 completed----")
            if l==int(len(content1)/4*3):
                print("----%75 completed----") 
            if l==int(len(content1)):
                print("----%100 completed----") 

    f.close()

                        #  ['FT', 'TOPO_DOM', '', '  6155..6306']


                    
                    
            
    print("---------getuniprot.py ended-----------\n\n")

