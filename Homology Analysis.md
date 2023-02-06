To find the homology between the retrieved sequences:

First a blast database is generated using these sequences
```
makeblastdb -in all_npxy.fasta -dbtype prot -parse_seqids -out npxy_blast_database.txt -title "my100proteins"
```
Than the the sequences are compared within this database to detect the homology. 
```
blastp -db npxy_blast_database.txt -num_alignments 0 -num_descriptions 2 -evalue 1e-20 -query all_npxy.fasta -out npxy_out.txt
```
The output file is quite unstructured and has to be parsed before analysis. The details of the parameters of the output can be found [here](https://web.cas.org/help/BLAST/topics/srch_sta.htm)

---

To recieve a structured data
```
blastp -db npxy_blast_database.txt -query all_npxy.fasta -out npxy_out.txt -outfmt "6 qseqid sacc qlen slen evalue bitscore pident qcovs"

sed -i '1s/^/qseqid sacc qlen slen evalue bitscore pident qcovs\n/' npxy_out.txt
```
parameters of this command can be found [here](https://www.biostars.org/p/88944/)

The generated file can be used to draw network graph

A cutoff for homology should be determined among the parameters

---


    