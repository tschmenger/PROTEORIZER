# PROTEORIZER
Proteorizer is an explorative tool that takes variants from laboratory or clinical settings and contextualizes the variants based on prior information from the protein of interest and similar proteins according to where these functional positions are located in the 3D structure of the protein of interest.
<img src="https://github.com/tschmenger/PROTEORIZER/blob/4e16b172e6d920daf4007e048936df43805d1d30/GraphicalAbstract.svg?sanitize=true">

# Please Cite
When using insights generated with Proteorizer please cite **PUBLICATION**.

# Introduction
This tool is a spritual successor to the work published in 2022 ```Schmenger et al. "Never-homozygous genetic variants in healthy populations are potential recessive disease candidates" (https://doi.org/10.1038/s41525-022-00322-z)```.

The idea of **Proteorizer** is to use information available for the protein of interest (your submission) and similar proteins (retrieved via Orthofinder or from your submitted custom alignment). Information of similar proteins is mapped back to the protein of interest, which is used to define functional clusters by measuring intramolecular distances.

Idea: If you know the function and effects of residues and mutations physically close to the position of interest, you might use this knowledge to better understand the effects of the variant in question.

We offer two different methods of how these clusters are defined ([Random Walk](https://igraph.org/r/doc/cluster_walktrap.html) and [Hierarchical Clustering](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/hclust)).

Several features were used to perform bayesian integration (see the paper [above](https://doi.org/10.1038/s41525-022-00322-z)). While **not being the main focus** of this approach we also provide scores based on the selected method. **Scores above 15** should be highlighting interesting cases.

The final results include the **functional and mapped-back information**, a **graphical representation of clusters** as well as an **annotated alignment**.


# How to use
## R-Shiny Application

### Preparing a custom alignment
The following steps can be used to create an alignment simply using blastp and clustal omega. <br>
**Note:** The bigger the alignment, the more time your request will consume. It is recommended to use not more than ~ 100 aligned sequences.

#### Step 1
Download the fasta sequence of your protein of interest. For RHOA you could do this via Uniprot, like [this](https://rest.uniprot.org/uniprotkb/P61586.fasta). <br>
**Note**: You can easily build this url using **https://rest.uniprot.org/uniprotkb/** + uniprotID +**.fasta**

#### Step 2
Use the downloaded fasta sequence to perform a blast search for similar sequences [blastp](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins). For more information on how to use BLAST please see [Blast Help](https://blast.ncbi.nlm.nih.gov/doc/blast-help/). **Make sure to select "Swissprot" as a database, or otherwise make sure that the retained accessions will be UniprotIDs.**

#### Step 3
Select the sequences you prefer and download them. <br>
**Note**: Make sure to download the complete sequences.
  
#### Step 4
Add the protein of interests fasta manually to the top of the just downloaded file, if it isn't present already. Change the formatting to roughly mimic the formatting of the remaining entries.
  
#### Step 5
Perform a multiple sequence alignment using [Clustal Omega](https://www.ebi.ac.uk/Tools/msa/clustalo/). Make sure to select **Protein**. Download & save the alignment file for usage with this script. Input the sequences via copy & paste or upload a file.

Make sure you download the complete MSA (including the clustal version, followed by 2 empty lines, followed by the MSA).

### Submitting your request
Please format the request following this simple rule:
**Identifier/Mutations**
- Identifier could be UniProtID or Gene Name. 
- Please separate multiple mutations using "," (comma).
Examples: 
- **P61586/Y34C,E40K** 
- **RHOA/Y34C**

## Command Line
### Requirements for command line usage
#### Python

#### R







# External (public) sources
