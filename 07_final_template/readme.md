## Final Project Template (currently not finished)

### 1. Get Ready

Install R on the computer you are working with if you have not already. Install any necessary packages as shown in 
the previous tutorials.

### 2. Download data

You will want to create a working directory (on your desktop, on a flash drive, all that matters is that you know what it is
called and where it is). In this working directory, create a new directory called: "raw_data" This directory will be where
you place all the fastq files for your project. 

Download #1: 

Download #2: 

### 3. Load packages one at a time to make sure they work.

Make sure after each "library" command that you do not get any errors otherwise you will run into problems later.

    library(dada2)
    
    library(phangorn)
    
    library(DECIPHER)
    
    library(phyloseq)
    
    library(ggplot2)
    
    library(reshape)
    
    library(RColorBrewer)
    
    library(plyr)
    
    library(picante)

### 4. Set up analysis

Set your working directory to that which you created. 

Use list.files() to check that you have the directory raw_data in your current path as well as the other files you downloaded
and put here. If not, your path is incorrect.

    # set the path to the location of the sequencing reads
    path <- "raw_data"
    
    list.files(path)
    
The files should have been listed with the list.files(path) command. If not, the directories are incorrect.

    # sort the order of the forward and reverse reads
    fnFs <- sort(list.files(path, pattern="_R1_001_trimmed.fastq"))
    
    fnRs <- sort(list.files(path, pattern="_R2_001_trimmed.fastq"))

    # extract sample names from files
    sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
    
    sample.names

You should be able to see the sample names from the second command above.

    # Specify the full path to the fnFs and fnRs
    fnFs <- file.path(path, fnFs)
    
    fnRs <- file.path(path, fnRs)
    
### 5. Error profiling

    # Look at a few sequences and check out their error profiles
    plotQualityProfile(fnFs[1:3])
    
    plotQualityProfile(fnRs[1:3])
    
    # file paths for putting the filtered reads
    filt_path <- file.path(path, "filtered")
    
    filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
    
    filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

    # filter and trim the samples
    out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(140,140),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=FALSE)
              
    head(out)
    
    # learn the error rates for the sequencing
    errF <- learnErrors(filtFs, multithread=FALSE)
    
    errR <- learnErrors(filtRs, multithread=FALSE)

    # look at plots of errors
    plotErrors(errF, nominalQ=TRUE)
    
    plotErrors(errR, nominalQ=TRUE)

### 6. Dereplicate and call sequence variants

    # dereplicate all reads
    derepFs <- derepFastq(filtFs, verbose=TRUE)
    
    derepRs <- derepFastq(filtRs, verbose=TRUE)
    
    # rename the dereplicated reads files
    names(derepFs) <- sample.names
    
    names(derepRs) <- sample.names

    # run the main file to call all of the unique sequences
    dadaFs <- dada(derepFs, err=errF, pool=T,multithread=TRUE)
    
    dadaRs <- dada(derepRs, err=errR, pool=T,multithread=TRUE)

Save progress:

    save.image(file="microbe_workflow1.Rdata")

### 7. Merge forward/reverse reads, and remove chimeras

    # merge the forward and reverse sequences
    mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
    
    # make a table of all sequences
    seqtab <- makeSequenceTable(mergers)
    
    # Inspect distribution of sequence lengths
    table(nchar(getSequences(seqtab)))
    
    # keep all mergers with length near the mode (253)
    seqtab <- seqtab[,nchar(colnames(seqtab)) %in% seq(252,254)]
    
    # remove chimeras
    seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
    
    # proportion of sequences not chimeric
    sum(seqtab.nochim)/sum(seqtab)

### 8. Summarize filtering

    # summarize the filtering
    getN <- function(x) sum(getUniques(x))
    
    track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
    
    colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
    
    rownames(track) <- sample.names
    
    track

### 9. Assign taxonomy

Here, we will use the GreenGenes database formatted for DADA2 to classify the 16S sequences we have here. Note that the
classification is only accurate to the family level, and any inferences about genus or species may or may not be accurate.

    taxa <- assignTaxonomy(seqtab.nochim, "gg_13_8_train_set_97.fa.gz", multithread=TRUE)
    
    unname(head(taxa))

Save progress:

    save.image(file="microbe_workflow2.Rdata")

### 10. Create phylogenetic tree of microbial data

    seqs <- getSequences(seqtab.nochim)
    
    names(seqs) <- seqs
    
    alignment <- AlignSeqs(DNAStringSet(seqs))
    
    phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
    
    dm <- dist.ml(phang.align)
    
    treeNJ <- NJ(dm)
    
    plot(treeNJ, show.tip.label=F)

Save progress:

    save.image(file="microbe_workflow3.Rdata")

### 11. Make a phyloseq object with the sample data, phylogeny, and sequence variant table

    # Make a data.frame holding the sample data
    samples.out <- rownames(seqtab.nochim)
    
You will have to make a custom data frame based on the sequences you have. For each characteristic of your data you are 
interested in, you should make a vector for each item first, then combine them all into a dataframe. An example is shown here 
with 10 samples, the sample IDs, and three additional traits. You should modify these, rename the traits,
add to or reduce vector size lengths as necessary given your sample size and the traits of interest
    
    sample.id <- sample.names
    
    trait1 <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j")
    
    trait2 <- c(1, 3, 1, 1, 5, 3, 5, 4, 4, 5)
    
    trait3 <- c("a", "a", "d", "d", "f", "f", "g", "g", "j", "j")
    
    micro.df <- data.frame(Sample.ID=sample.id, trait_name1=trait1, trait_name2=trait2, trait_name3=trait3)

    rownames(micro.df) <- samples.out


    # construct a phyloseq object
    ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=F),
    		sample_data(micro.df),
    		tax_table(taxa),
    		phy_tree(treeNJ))
        
    ps
    
    # remove cyanobacteria
    ps <- subset_taxa(ps, Phylum != "p__Cyanobacteria")

Save progress:

    save.image(file="microbe_workflow4.Rdata")

### 12. Composition

Look at a summary of the numbers of unique sequence variants:
    
    summary(sample_sums(ps))

Find the unique phyla, classes, and families in the dataset:

    sort(get_taxa_unique(ps, "Phylum"))
    
    sort(get_taxa_unique(ps, "Class"))
    
    sort(get_taxa_unique(ps, "Family"))

Plot the composition (by phylum)

    plot_bar(ps, x="Sample.ID", fill="Phylum") + facet_wrap(~Location, scales="free_x")

### 13. Alpha diversity

Make a table of many types of alpha diversity:

    alpha_table <- estimate_richness(ps)
    
    alpha_table

Let's make a simple table that includes 4 measures of alpha diversity. We'll save them to the object named 'alpha.'
	
    alpha <- cbind(estimate_richness(ps)[,c(1,6,7)], pd((otu_table(ps)@.Data), ps@phy_tree, include.root=F)$PD)
  
    colnames(alpha) <- c("Observed_SVs", "Shannon", "Simpson", "PD")
  
    alpha

Remember that sampling depth can influence estimates of diversity, especially if they are incompletely sampled communities.
Let's plot some rarefaction curves to take one look at this factor:

    rarecurve(otu_table(ps), step=5, cex=0.5, label=F, xlim=c(0,20000))

Let's plot a couple types of alpha diversity by trait3 (you can use the name of any of your traits here (if they are 
categorical):

    plot_richness(ps, x="trait3", measures=c("Shannon", "Simpson"))
     
We can also plot the relationships among each of the estimates of alpha diversity. Here, Observed_SVs = number of observed
sequence variants and PD = phylogenetic diversity.

    plot(alpha, pch=19, cex=1)

### 14. Beta diversity

Now we'll measure beta diversity among communities in a few different ways, colored by trait3 (again, any of your traits will
work. 

First, we'll make a plot of the Bray-Curtis distance:
	
    ord.nmds.bray <- ordinate(ps, method="NMDS", distance="bray")

    plot_ordination(ps, ord.nmds.bray, color="trait3", shape="trait3", title="Bray NMDS") + geom_point(size=3)

And unweighted Unifrac distance:

    ord.unweighted.unifrac <- ordinate(ps, "PCoA", "unifrac", weighted=F)

    plot_ordination(ps, ord.unweighted.unifrac, color="trait3", shape="trait3",title="Unweighted Unifrac PCoA") + geom_point(size=3)
	
And weighted Unifrac distance:

    ord.weighted.unifrac <- ordinate(ps, "PCoA", "unifrac", weighted=T)

    plot_ordination(ps, ord.weighted.unifrac, color="trait3", shape="trait3",title="Weighted Unifrac PCoA") + geom_point(size=3)

Save progress:

    save.image(file="microbe_workflow5.Rdata")

### 15. Stats

Use a relevant statistical test that was covered in class (or if you want, one we haven't covered) to test your hypotheses.
This will be different for each dataset and question, and will likely involve a little bit of trouble shooting. Feel free
to ask questions or search the internet for solutions. I expect this may take up more time than many of the other parts
of the project.
