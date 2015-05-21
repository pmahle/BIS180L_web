---
layout: lab
hidden: true
title: 'Metagenomics'
output: html_document
tags:
- Linux
- Illumina
- Metagenomics
---

## Speaking notes
metagenomics  
16s rRNA  
Rice root  
OTUs  
Alpha vs beta diversity  

## Kristen's testing stuff

+ Right click "Extract to "QIIME-1.9.0-amd64.vdi\" (unless Tim unzipped it)
+ When unzipped it's 10.7 GB - 15 minutes transfer via USB3 (much faster!!)
	+ When zipped 3:35 transfer and to unzip on the drive it took 33 minutes
+ The student's current USB is 10-14.5GB

## QQ for Julin
- Nipponbare early vs late? is this young and older plants? or time of year? or time of day?

____
# Assignment

Pull your repository for the assignment template. 

Knit the file and submit the .Rmd and .html when you are done.  Open an issue indicating that the assignment is ready to be graded.

# Background

Metagenomics is a rapidly expanding field that has the power to explain microbial communities with a very high resolution by leveraging next generation sequencing. There are applications in the clinic, in ecological environments, food safety, and others. By definition, metagenomics is the study of a collection of genetic material (genomes) from a mixed community of organisms typically microbial.  

Today we will walk through a common metagenomics workflow using QIIME (pronounced "chime") by completing the following:

1. Determine the various microbial communities in our samples
2. Calculate the diversity within our sample (alpha diversity)
3. Calculate the diversity between different sample types (beta diversity)

*Acknowledgement must be paid to Professor Scott Dawson for sharing his original metagenomics lab that we have adapted for this class*

## Getting Started with QIIME
Quantitative Insights Into Microbial Ecology or QIIME is an open-source bioinformatics pipeline for performing microbiome analysis from raw DNA sequencing data. It has been cited by over 2,500 peer-reviewed journals since its [publication](http://www.nature.com/nmeth/journal/v7/n5/full/nmeth.f.303.html) in 2010.  

QIIME requires many dependencies which can make installing it a bit of headache. However, the developers of QIIME have made a standalone Virtual Box with a complete install. Even though, we know you're install pros we will be working with this more convenient setup. The QIIME Virtual Box has been downloaded to the desktop in SLB2020

+ Plug in your USB to a lab computer
+ Transfer the folder titled QIIME to your USB by dragging and dropping (or copy and pasting).
+ This will take 15 minutes or so. 
+ The **password** for the QIIME Virtual Box is qiime
+ QIIME is very RAM heavy. Therefore, some of the steps you'll complete today will take a couple of minutes to run. Please be patient.


## Background for our Data
Today, we will be working with the samples collected from the rhizosphere of rice plants. The rhizosphere is an area of soil near the plant roots that contains both bacteria and other microbes associated with roots as well as secretions from the roots themselves. See diagram below from [Phillppot et al., *Nature*, 2013](http://www.nature.com/nrmicro/journal/v11/n11/full/nrmicro3109.html).
![plot of rhizosphere]({{ site.baseurl }}/figure/metagenomics_lab-1-rhizosphere.jpg) 

Samples were collected and sequenced XXX TODO. We will be working with the sequencing results in `RiceSeqs.fna` and sample information in `RiceMappingFile.txt`. These are already downloaded to the QIIME Virtual Box in `~/Desktop/Data`.
16S rRNA profiling

From Julin: 454 sequence data (remember them). 16S sequences from samples taken within the root, at the root surface, near the root, and then in soil further away. There should be a clear spatial signal.

The data we will be using for this lab has already been installed on the desktop of the QIIME virtual box. The barcodes have been removed. TODO add in more info about the data set.

## Explore and Quality Control Data
Often times, as bioinformaticians, we will receive data sets with little background. Sometimes the first step is to explore the raw data that we will be working with. This can help us spot inconsistencies or logical fallacies down the line when working with more automated pipelines.  

Open RiceMappingFile.txt with `less` to view more information about the data you are working with. This file contains information about each sample including the cultivar, treatment, and number of technical replicates. It also includes the barcodes used to identify each sample during multiplexing which should be in a 1:1 ratio. Let's use the barcodes to determine if we have an even number of reads per sample type.  

In the RiceSeqs.fna file, barcodes for each sequence are indicated in the header with `new_bc=`. These barcodes are also mapped to the sample information in RiceMapping.txt. Try to determine the number of sequences present for each barcode. This can be accomplished using just Linux/Unix commands. I'll start by giving you the tools, so you can try to piece together the command on your own. 

**Helpful Commands (in no particular order):** `cut`, `grep`, `head`, `sort`, `uniq` and good 'ol `|` to chain the commands together.

If you get stuck, highlight the hidden text underneath this sentence for one potential solution.  
<font color="white" face="menlo">
grep ">" RiceSeqs.fna | cut -d " " -f 4 | sort | uniq -c
</font>

**Exercise 1**
Using information in the RiceMappingFile.txt and RiceSeqs.fna answer the following questions. Are the number of sequences for each sample type approximately the same or are there any outliers? If so, which samples do they belong to? Could a different number of sequences per sample affect future analysis? Explain your reasoning.

Exercise 1 KEY
Yes, outliers
IM1& IM2 are highest - IR50, 1mm_soil, tech rep2
MB1 & MB2 have the lowest - root, M04 cultivar, root surface

Now that we've poked around in our raw data a little. Let's carry on with analyzing the microbes present in our samples.

## Classify Various Microbiome Sequences into OTUs
Operational taxonomical units (OTUs) are used to describe the various microbial species in a sample. OTUs are defined as a cluster of reads with 97% 16S rRNA sequence identity. We will use QIIME to classify OTUs into an OTU table.

```bash
cd ~/Desktop 
pick_de_novo_otus.py -i Data/RiceSeqs.fna -o otus_demo
# use ctrl-c to stop this process, see note below
```
`-i` designates the input file with our DNA sequences  
`-o` designates the folder for our output  
The syntax to designate input and output files is consistent across QIIME commands.
**Note:** this command takes ~45 minutes to run and uses all RAM for the Virtual Box basically freezing the system until it's complete so we have ran this for you and added the data in `~/Desktop/otus`.

Let's view the statistics of the OTU table (a binary file) and save the output.

```bash
biom summarize-table -i otus/otu_table.biom > otus/otu_class.txt
```
**Exercise 2**
From the OTU summary, look at how many OTUs correspond to each sample ("counts/sample detail"). Do technical replicates agree with one another? At this stage, can you observe any trends in number of OTUs based on cultivar or sample location?  
*Note:* The OTUs actually match 1:1 with the number of sequences per sample ID/barcode. This is an artifact of creating the smallest demo data set to use on the Virtual Boxes. The OTUs are still representative of the microbial diversity in that sample though.

Exercise 2 KEY
Scott's original question  
What is the minimum/maximum of OTUs in all samples? Are there significant differences between the samples? If so, why do you think?  
Min = IE1 289  
Max = IM2 4880  
36155 total OTUs  
OTU assignment is mutually exclusive (since such a high percent identity)  

As we learned last week, we can rely on the human eye to help pick out patterns based on color. We are going to make a heat map of the OTUs per sample. The OTU table is visualized as a heat map where each row corresponds to an OTU and each column corresponds to a sample. The higher the relative abundance of an OTU in a sample, the more intense the color at the corresponding position in the heat map. OTUs are clustered by [UPGMA hierarchical clustering](http://en.wikipedia.org/wiki/UPGMA). QIIME indicates the biological classification by prefixing the level for example "p\_" indicates phylum and  "g\_" indicates genus. For a refresher on biological classification, view this [helpful wiki page](http://en.wikipedia.org/wiki/Bacterial_taxonomy).  

```bash
make_otu_heatmap.py -i otus/otu_table.biom -o otus/OTU_Heatmap.pdf
```

**Exercise 3**
Although, the resolution of the y-axis makes it different to read each OTU, it is still a valuable preliminary visualization. What types of information can you gain from this heat map? Are there any trends present at this stage with respect to the various samples?

Now we'd like to visualize our data with a little higher resolution and summarize the communities by their taxonomic composition.

```bash
summarize_taxa_through_plots.py -i otus/otu_table.biom -o otus/wf_taxa_summary -m Data/RiceMappingFile.txt
```
`-m` provides the path to the mapping file with sample meta data

Summarize OTU by Category (optional, pass -c); Summarize Taxonomy; and Plot Taxonomy Summary

From the bar charts, which groups are the predominant phyla in the different samples? Are there any predominant groups unique to particular samples? Do you have any explanations for any of the observed differences between the different samples?


5. Compute alpha diversity-generate rarefaction curves

```bash
echo “alpha_diversity:metrics shannon,PD_whole_tree,chao1,observed_species” > otus/alpha_params.txt
alpha_rarefaction.py –i otus/otu_table.biom –m Data/RiceMappingFile.txt –o otus/wf_arare -p otus/alpha_params.txt –t otus/rep_set.tre
```
First, we make the file that tells contains the parameters to make the alpha rarefaction. Then we calculate it.
Lots of helpful info in this plot

6. Compute beta diversity and generate PCOA plots and UPGMA trees

PCOA plots:
```bash
beta_diversity_through_plots.py –i otus/otu_table.biom –m RiceMappingFile.txt –o wf_bdiv_even289 -t otus/rep_set.tre –e 289
```
You may see a negative Eigenvalue error, but the negative values are 100x smaller than the positive values so we can ignore the warning.
QQ - what is 289 doing?
QQ - is there a 2D output do we need a different option?

UPGMA trees
```bash
￼upgma_cluster.py –i wf_bdiv_even289/unweighted_unifrac_dm.txt –o
￼unweighted_upgma.tre
￼upgma_cluster.py –i wf_bdiv_even289/weighted_unifrac_dm.txt –o
￼weighted_upgma.tre
```
The files saved here are the distances used to form a tree of the samples
Upload them to http://iubio.bio.indiana.edu/treeapp/treeprint-form.html
"Phenogram" will be the most useful display







Use MicrobeWiki (or other online sites to learn more about these microbes) o http://microbewiki.kenyon.edu/index.php/MicrobeWiki

_____

## Annotate the differentially expressed genes.

As part of [his B_rapa annotation paper](http://www.g3journal.org/content/4/11/2065.long) a postdoc in my lab, Upendra Devisetty, did a quick annotation of the _B. rapa_ genes by BLASTing each gene to the NCBI non-redundant database and taking the gene description of the best match (above a threshold). You can download the description with this command:

    wget http://www.g3journal.org/content/suppl/2014/08/12/g3.114.012526.DC1/FileS9.txt
  
Place the file into your `Brapa_reference` directory.
  
Now open Rstudio, set the working directory to your `Diff_Exp` directory, and proceed.
  
__Exercise 1__

__a.__ Use `merge()` to add gene descriptions for the genes found to be regulated by the DP treatment.  Output a table of the top 10 genes that includes the output from edgeR and the descriptions.  __Important: Pay attention to the "sort="" argument to `merge()`.  Should it be TRUE or FALSE?

__b.__ Repeat this for  genes with a genotype x trt interaction.



By looking at the annotated interaction gene list we can see that many of the identified genes code for proteins that modify the plant cell wall (I wouldn't expect you to know this unless you are a plant biologist).  This might relate to the different properties of the two cultivars, IMB211 and R500, and their different responses to the treatment.  Depending on our interests and knowledge of the system we might at this point choose follow up study on specific genes.

## Test for enrichment of functional classes of genes.

A casual glance indicated that there might be an enrichment for cell wall related genes in  the gt X trt DE gene list.  We can test this more rigorously by asking if there is statistical over-representation of particular [Gene Ontology (GO)](http://geneontology.org/) terms in this set of genes.  GO terms provide a precise, defined language to describe gene function.

To test for test for enrichment of GO terms we essentially ask the question of whether the prevalence of a term in our gene set of interest is higher than its prevalence in the rest of the genome (aka the universe).  For example if 20% of the differentially expressed genes have the GO term "Cell Wall" but only 10% of the not-differentially expressed genes have the term Cell Wall that the term "Cell Wall" might be over-represented, or enriched, in our differentially expressed genes.  In principle this could be tested using a [Chi-squared test](http://www.biostathandbook.com/chiind.html) or [Fisher's exact test](http://www.biostathandbook.com/fishers.html) for contingency tables.  In practice we will use [GOseq](http://www.bioconductor.org/packages/release/bioc/html/goseq.html) that [makes an adjustement for gene-length bias](http://genomebiology.com/2010/11/2/r14) since it is easier to detect differential expression of longer genes.  

It is important in these analyses to define the "Universe" correctly.  The Universe of genes is all of the genes where you could have detected an expression difference.  It would not have been possible to detect an expression difference in genes that were not expressed at a detectable level in our experiment.  So the Universe in this case is all expressed genes, rather than all genes.

### Install GOseq


```r
source("http://bioconductor.org/biocLite.R")
biocLite("goseq")
```

### Download the GO annotation and gene length for B.rapa

GO term annotation for _B. rapa_ is also available from the [Devisetty et al paper](http://www.g3journal.org/content/4/11/2065.long).

cDNA gene lengths were estimated by the `featureCounts()` command used in Tuesday's lab.  I have saved them for you and you can download them as instructed below.  (You're welcome)

Download the files using the commands below (in Linux) and place the file in your `Brapa_reference` directory.

GO terms:

    wget http://www.g3journal.org/content/suppl/2014/08/12/g3.114.012526.DC1/FileS11.txt
    
Gene Length:

    wget http://jnmaloof.github.io/BIS180L_web/data/Brapa_CDS_lengths.txt
  
### Format data for GOseq

We need to do a bit of data wrangling to get things in the correct format for GOSeq.  First we import the gene lengths and GO terms.  We also import a table of all expressed genes in the experiment (you could have gotten this from Tuesday's data)

Get the data

```r
library(goseq)
go.terms <- read.delim("../Brapa_reference/FileS11.txt",header=FALSE,as.is=TRUE)
head(go.terms)
names(go.terms) <- c("GeneID","GO")
summary(go.terms)

expressed.genes <- read.delim("internode_expressed_genes.txt",as.is=TRUE)
head(expressed.genes)
names(expressed.genes) <- "GeneID"

gene.lengths <- read.table("../Brapa_reference/Brapa_CDS_lengths.txt",as.is=TRUE)
head(gene.lengths)
summary(gene.lengths)

#we need to reduce the gene.length to only contain entries for those genes in our expressed.genes set.  We also need this as a vector
gene.lengths.vector <- gene.lengths$Length[gene.lengths$GeneID %in% expressed.genes$GeneID]
names(gene.lengths.vector) <- gene.lengths$GeneID[gene.lengths$GeneID %in% expressed.genes$GeneID]
head(gene.lengths.vector)

#Do the reverse to make sure everyting matches up (it seems that we don't have length info for some genes?)
expressed.genes.match <- expressed.genes[expressed.genes$GeneID %in% names(gene.lengths.vector),]
```

Format go.terms for goseq.  We want them in list format, and we need to separate the terms into separate elements.

```r
go.list <- strsplit(go.terms$GO,split=",")
names(go.list) <- go.terms$GeneID
head(go.list)
```

Format gene expression data for goseq.  We need a vector for each gene with 1 indicating differential expression and 0 indicating no differential expression.

```r
DE.interaction <- expressed.genes.match %in% rownames(DEgene.interaction) 
    #for each gene in expressed gene, return FALSE if it is not in DEgene.trt and TRUE if it is.
names(DE.interaction) <- expressed.genes.match
head(DE.interaction)
DE.trt <- as.numeric(DE.interaction) #convert to 0s and 1s
head(DE.interaction)
sum(DE.interaction) # number of DE genes
```

### Calculate over-representation

Now we can look for GO enrichment

```r
#determines if there is bias due to gene length.  The plot shows the relationship.
nullp.result <- nullp(DEgenes = DE.interaction,bias.data = gene.lengths.vector)

#calculate p-values for each GO term
rownames(nullp.result) <- names(gene.lengths.vector) #because of a bug in nullp()
GO.out <- goseq(pwf = nullp.result,gene2cat = go.list,test.cats=("GO:BP"))

#list over-represented GO terms (p < 0.05)
GO.out[GO.out$over_represented_pvalue < 0.05,]
```

### GO visualization

Looking through a long list can be tough.  There is a nice visualizer called [REVIGO](http://revigo.irb.hr/).  To use it we need to cut and paste the column with the GO term and the one with the p-value.  Use the command below to print these to columns to the console:


```r
print(GO.out[GO.out$over_represented_pvalue < 0.05,1:2],row.names=FALSE)
```

Cut and paste the GO terms and p-values into [REVIGO](http://revigo.irb.hr/).  You can use the default settings.  There are three types of GO terms:

* Biological Process (BP)
* Cellular Compartment (CC)
* Molecular Functions (MF)

Generally I find the BP terms most helpful but you can look at each type by clicking in the tabs.

REVIGO has three types of visualizers.  The "TreeMap"  Can be particular nice because it groups the GO terms into a hierarchy.

__Exercise 2__:  

__a:__ In REVIGO display a "TreeMap" of the BP GO terms.  Was our hypothesis that cell wall genes are enriched in the genotype X treatment gene set correct?  You DO NOT need to include the treemap in your answer.

__b:__ Display a "TreeMap" of the CC GO terms.  There are four general categories shown, some with sub-categories.  What are the two general categories with the largest number of sub categories?  How might these general categories relate to differences in plant growth?  You DO NOT need to include the treemap in your answer.

## Promoter motif enrichment

To help us understand what causes these genes to be differentially expressed it would be helpful to know if they have any common transcription factor binding motifs in their promoters.

First we must get the sequence of the promoters.  For ease we will define the promoters as the 1000bp upstream of the start of the gene.  In the interest of time I have pre-computed this for you, but I show you how to do this below in case you need to do it in the future.

### Get gene "promoters".

__You do not need to run this__

To extract the promoters I used Mike Covington's [extract-utr script](https://github.com/mfcovington/extract-utr) and downloaded the CDS as file S4 from [Devisetty et al](http://www.g3journal.org/content/4/11/2065.long)

The command I used was
```
extract-utr.pl --gff_file=Brapa_gene_v1.5.gff \
--genome_fa_file=BrapaV1.5_chrom_only.fa  \
--cds_fa_file=Brassica_rapa_final_CDS.fa  \
--fiveprime --utr_length=1000 --gene_length=0 \
--output_fa_file=BrapaV1.5_1000bp_upstream.fa
```

### Gather data for motif enrichment

First lets gather the data that we need.  you do need to run the following.

Download the promoters and place them in your `Brapa_reference` directory.  You can download them with

    wget http://jnmaloof.github.io/BIS180L_web/data/BrapaV1.5_1000bp_upstream.fa.gz
  
Siobhan Brady has compiled a file of plant transcription factor binding motifs.  You can download those from

    wget http://jnmaloof.github.io/BIS180L_web/data/element_name_and_motif_IUPACsupp.txt
  
Place these in your `Brapa_reference` directory

Load the promoter sequences

```r
library(Biostrings) #R package for handling DNA and protein data
promoters <- readDNAStringSet("../Brapa_reference/BrapaV1.5_1000bp_upstream.fa.gz")

#convert "N" to "-" in promoters.  otherwise motifs will match strings of "N"s
promoters <- DNAStringSet(gsub("N","-",promoters))

promoters
```

Load the motifs and convert to a good format for R

```r
motifs <- read.delim("../Brapa_reference/element_name_and_motif_IUPACsupp.txt",header=FALSE,as.is=TRUE)
head(motifs)
motifsV <- as.character(motifs[,2])
names(motifsV) <- motifs[,1]
motifsSS <- DNAStringSet(motifsV)
motifsSS
```

Next we need to subset the promoters into those in our DE genes and those in the "Universe"

```r
#get names to match...there are are few names in the DEgene list not in the promoter set
DEgene.interaction.match <- row.names(DEgene.interaction)[row.names(DEgene.interaction) %in% names(promoters)]

#subset promoter files
universe.promoters <- promoters[expressed.genes.match]
target.promoters <- promoters[DEgene.interaction.match]
```

### Look for over-represented motifs

I have written a function to wrap up the required code.  Paste this into R to create the function

```r
#create a function to summarize the results and test for significance
motifEnrichment <- function(target.promoters,universe.promoters,all.counts=F,motifs=motifsSS) {
  
  #use vcountPDict to count the occurences of each motif in each promoter
  target.counts <- vcountPDict(motifs,target.promoters,fixed=F) + 
    vcountPDict(motifsSS,reverseComplement(target.promoters),fixed=F)
  universe.counts <- vcountPDict(motifs,universe.promoters,fixed=F) + 
    vcountPDict(motifsSS,reverseComplement(universe.promoters),fixed=F)
  
  if (all.counts) { 
    #count all occurences of a motif instead of the number of promoters that it occurs in
    target.counts.sum <- apply(target.counts,1,sum)
    universe.counts.sum <- apply(universe.counts,1,sum)
  } else {
    target.counts.sum <- apply(ifelse(target.counts > 0,1,0),1,sum)
    universe.counts.sum <- apply(ifelse(universe.counts > 0 , 1, 0),1,sum)
  }
  n.motifs <- length(target.counts.sum)
  results <- vector(mode="numeric",length=n.motifs)
  for (i in 1:n.motifs) {
    if (all.counts) { #the contigency tables are different depending on whether we are looking at promoters or overall occurences
      #test if ratio of occurences to promoters is the same in the target and the universe
      m <- matrix(c(
        target.counts.sum[i],                       #number of occurences within target
        dim(target.counts)[2],                      #number of promoters in target
        universe.counts.sum[i],                  #number of occurences within universe
        dim(universe.counts)[2]                  #number of promoters in universe
      ),ncol=2)
    } else { #looking at promoters with and without hits
      m <- matrix(c(
        target.counts.sum[i],                        #number of promoters in target with hit
        dim(target.counts)[2]-target.counts.sum[i],            #number of promoters in target with no hit
        universe.counts.sum[i],                   #number of promoters in universe with hit
        dim(universe.counts)[2]-universe.counts.sum[i]   #number of promoters in universe with no hit
      ),ncol=2)
    } #else
    results[i] <- fisher.test(m,alternative="greater")$p.value
  } #for loop
  results.table <- data.frame(
    motif=names(motifs),
    universe.percent = round(universe.counts.sum/dim(universe.counts)[2],3)*100,
    target.percent = round(target.counts.sum/dim(target.counts)[2],3)*100,
    p.value =  results)
  results.table <- results.table[order(results.table$p.value),]
  results.table
}
```

Now with the function entered we can do the enrichment

```r
motif.results <- motifEnrichment(target.promoters,universe.promoters)
head(motif.results)
```
The resulting table gives the p-value for enrichment for each motif, as well as the %of promoters in the universe and in our target gene set that have the motif.


__Exercise 3__ 

__a.__ How many motifs are enriched at P < 0.05?  
__b.__ What is the identity of the most significantly over-enriched promoter?  
__c.__ What percentage of genes in the "Universe" have this motif?  What percentage in our target set?  
__d.__ You can find information on the motifs [here](http://arabidopsis.med.ohio-state.edu/AtcisDB/bindingsites.html).  Do you think that the most enriched motif represents a biologically meaningful result?  Discuss why or why not.









