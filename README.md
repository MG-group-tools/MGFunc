# MGFunc - Metagenomics Functional and Taxonomical annotation Pipeline

Authors:
Asli Ismihan Ozen (asli@cbs.dtu.dk)
Kosai Al-Nakeeb (kosai@cbs.dtu.dk)
Thomas Sicheritz-Ponten (thomas@cbs.dtu.dk)

MGFunc is an automated pipeline for annotating metagenomic gene catalogues using remote homology based functional assessment and phylogenetic analysis for function and taxonomy. It provides a workflow solution on protein level to process large metagenomics data resulting from high-throughput sequencing experiments. MGFunc offers the researchers a customizable tool to identify, cluster and analyze the protein space. The main input to the pipeline is a set of unknown protein sequences and the outcome is a list of annotated ortholog clusters and novel ortholog clusters that are presumably functional homologs.

# Input
A Metagenomics gene catalog colelcted from one or several samples(preferably all samples to be analyzed)

# Methods
1. UBLAST 
2. BLASTDECIDE
3. CLUSTERING
4. FILTERING
5. PARSING DATABASE INFO
6. EXTRACTION
7. SPLIT & CHANGE HEADERS
8. ALIGNMENT
9. TREES

# Output
Gene clusters with functional agreement, gene function annotation for each gene with GO terms, E.C. numbers if available, KEGG pathway links  and taxonomical annotation throuugh closest relative.



## METHODS 
### 1. UBLAST 
Gene sequence search in the protein datbaases. Currently in Uniprot Database. 
Internal homology detection within the gene catalog, to use in the clustering section.

### 2. BLASTDECIDE
Selection of all the significant hits for a given gene. 

### 3. CLUSTERING
Clustering genes based on reciprocal hits principle. Genes that have reciprocal hits are connected in a directional graph and the connected components of the graph are considered as a
cluster. 

### 4. FILTERING
Choose clusters to work with based on size or how many different samples they inherit, or specific samples they inherit.

### 5. PARSING DATABASE INFO
Getting all the database sequence hits for the genes in each cluster and collectiong database metadata about the cluster genes. And checks the functional consistency within a group. 

### 6. EXTRACTION
Extract sequences for each cluster from both gene catalog and corresponding database hits, saving in a fasta format where headers represent the cluster ID and uniqe gene id. 

### 7. SPLIT & CHANGE HEADERS
Split sequences for different clusters from the same file and change the fasta seq headers to short names due to the further use in tree making. Saves a conversion table

### 8. ALIGNMENT
Muscle or mafft alignment, fast or slower options can be chosen. 

### 9. TREES
Making trees using PAUP. 


