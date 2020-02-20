## Microbiome Community Abundance and Diversity 20 February 2020

We will continue from the preliminary filtering of microbiome data we previously did in class. Here, you will
be plotting and interpreting information about the communities' composition and alpha and beta diversity measures.

Make sure you are in the same working directory as before with the data and the saved ".RData" file. 

### 1. Load packages one at a time to make sure they work

    library(dada2)
    
    library(phangorn)
    
    library(DECIPHER)
    
    library(phyloseq)
    
    library(ggplot2)
    
    library(reshape)
    
    library(RColorBrewer)
    
    library(plyr)
    
    library(picante)

### 2. Load your saved image of the work you did previously

    load("microbe_workflow1.Rdata")
    
    # check that the items loaded
    ls()

### 3. Create phylogenetic tree of microbial data

    seqs <- getSequences(seqtab.nochim)
    
    names(seqs) <- seqs
    
    alignment <- AlignSeqs(DNAStringSet(seqs))
    
    phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
    
    dm <- dist.ml(phang.align)
    
    treeNJ <- NJ(dm)
    
    plot(treeNJ, show.tip.label=F)

Save progress:

    save.image(file="microbe_workflow_day2.Rdata")

### 4. Make a phyloseq object with the sample data, phylogeny, and sequence variant table

    # Make a data.frame holding the sample data
    samples.out <- rownames(seqtab.nochim)
    samp.number <- sapply(strsplit(samples.out, "_"), `[`, 1)
    sample.id <- sapply(strsplit(samples.out, "_"), `[`, 2)
    genus.id <- sapply(strsplit(samples.out, "_"), `[`, 3)
    species.id <- sapply(strsplit(samples.out, "_"), `[`, 4)
    location <- sapply(strsplit(samples.out, "_"), `[`, 5)
    species.location <- paste(species.id, location, sep=".")
    micro.df <- data.frame(Sample.number=samp.number, Sample.ID=sample.id, Genus.ID=genus.id,
						Species.ID=species.id, Location=location, Species.Location=species.location)
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

    save.image(file="microbe_workflow_day2.Rdata")

### 5. Composition

Look at a summary of the numbers of unique sequence variants:
    
    summary(sample_sums(ps))

Find the unique phyla, classes, and families in the dataset:

    sort(get_taxa_unique(ps, "Phylum"))
    sort(get_taxa_unique(ps, "Class"))
    sort(get_taxa_unique(ps, "Family"))

Sample a random family:

	  sample(sort(get_taxa_unique(ps, "Family")), 1)

Plot the composition (by phylum)

    plot_bar(ps, x="Sample.ID", fill="Phylum") + facet_wrap(~Location, scales="free_x")

### 6. Alpha diversity

Here, we'll take a look at a couple ways you can investigate alpha diversity. There is a base function in the phyloseq 
package to look at many types of alpha diversity:

	estimate_richness(ps)

However, we'll just be looking at a few to keep things simple. We'll save them to the object named 'alpha.'
	
	alpha <- cbind(estimate_richness(ps)[,c(1,6,7)], pd((otu_table(ps)@.Data), ps@phy_tree, include.root=F)$PD)
  
	colnames(alpha) <- c("Observed_SVs", "Shannon", "Simpson", "PD")
  
	alpha

Remember that sampling depth can influence estimates of diversity, especially if they are incompletely sampled communities.
Let's plot some rarefaction curves to take one look at this factor:

	rarecurve(t(otu_table(ps)), step=5, cex=0.5, label=F, xlim=c(0,4000))

If we had the full datasets, the picture may look more completely sampled. Anyways, there are ways of using the full dataset
and also ways of rarefying prior to estimating alpha diversity and beta diversity. We will not get into those methods in these
activities.

Let's plot a couple types of alpha diversity:

	plot_richness(ps, x="Location", measures=c("Shannon", "Simpson"))
  
Now, let's reduce the sampling to be equal in each community and plot the diversity:
  
    rarefied <- rarefy_even_depth(ps, 200)
    
    plot_richness(rarefied, x="Location", measures=c("Shannon", "Simpson"))
     
We can also plot the relationships among each of the estimates of alpha diversity. Here, Observed_SVs = number of observed
sequence variants and PD = phylogenetic diversity.

	plot(alpha, pch=19, cex=1)

### 7. Beta diversity

Now we'll measure beta diversity among communities in a few different ways. 

First, we'll make a plot of the Bray-Curtis distance:
	
	ord.nmds.bray <- ordinate(ps, method="NMDS", distance="bray")

	plot_ordination(ps, ord.nmds.bray, color="Species.Location", shape="Species.Location", title="Bray NMDS") + geom_point(size=3)

And unweighted Unifrac distance:

	ord.unweighted.unifrac <- ordinate(ps, "PCoA", "unifrac", weighted=F)

	plot_ordination(ps, ord.unweighted.unifrac, color="Species.Location", shape="Species.Location",title="Unweighted Unifrac PCoA") + geom_point(size=3)
	
And weighted Unifrac distance:

	ord.weighted.unifrac <- ordinate(ps, "PCoA", "unifrac", weighted=T)

	plot_ordination(ps, ord.weighted.unifrac, color="Species.Location", shape="Species.Location",title="Weighted Unifrac PCoA") + geom_point(size=3)


Save progress:

    save.image(file="microbe_workflow_day2.Rdata")
