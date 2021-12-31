# EVOLUTIONARY CONSERVED MOTIFS ✨
##### _This tool is in development process_
#
#### Installation


>The packages used are;
>os,glob,matplotlib,scipy,pandas,numpy
#
```sh
git clone git@github.com:ZeynepKilinc/bbb-project.git
```
#### Run
```sh
python3 controller.py
```
The controller.py will automatically run other functions. This script currentlu needs manually entered motif directly inside the script _-this part will be fixed soon-_. 

## HOW THE PROCESS WORKS
#### ----------------------------------------------

#### Retrieving the proteins
- User provides a motif _-will be updated-_
- Motif will be searched in the ExPASY database, for Homo Sapiens. A protein list having this motif will be retrieved
- The proteins in the retrieved list will be searched in UniProt database.
- Proteins that are not having a transmembrane region or having the motif in their cytosolic domain will be eliminated.
- The remaining protein's sequence, database codes and position information of the motif will be retrieved and saved in the folder {motif}_proteins, as {motif}_proteins.tsv

>Some proteins may have the motif in the same/different cytosolic region for more than once. In such cases each occurence of the motif will be written seperately in the file in different lines. _see N-P-X-Y\_proteins.tsv_

#### Getting orthologs
- {motif}_orthologs folder will be created
- Using the provided OMA id's in the UniProt database, ortholog sequences of each protein will be retrieved.
- Retrieved sequence of the referenced protein in the previous step, will be added to the fasta file of its orthologs.
- Fasta files will be written in the created ortholog folder with the heading {uniProtId}_orthologs.fasta

#### Making the alignments
- {motif}_msa folder will be created
- MSA files will be created for each {uniProtId}_orthologs.fasta file, using ClustalO
- Aligned files will be stored in the msa folder with the heading {uniProtId}_orthologs_msa.fasta

#### Getting the conservations
