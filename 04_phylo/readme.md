## Some Basic Methods Manipulating Phylogenies 25 February 2020

Today we will use R to manipulate phylogenies. This may be useful if you wish to see if host phylogeny is related with
microbiome community composition. For example, you may wish to test if host phylogenetic distance is related to beta
diversity measurements. 

### 1. Set up working directory and download data

Make a directory on your desktop called 'phylogeny' and set that as your current working directory in RStudio. 

Download an example phylogeny from this link: [link](https://drive.google.com/file/d/1UKYfkm2TzM9B3lusF63sOeumjX8TGI6v/view).
Put this file in your working directory. This phylogeny is from whole-genome data of the bird _Certhia americana_. This is 
not really relevant to the exercise, but only if you were curious about the phylogeny you were working with.

### 2. Install necessary packages

    install.packages("ape")
    
    install.packages("phytools")

### 3. Load packages one at a time to make sure they work

    library(ape)
    
    library(phytools)

### 4. Working with a small fake phylogeny for initial practice.

First, we will create a fake phylogeny with 5 tips named A through E and call it "tree": 

    tree <- read.tree(text = "(((A,B),(C,D)),E);")

If we plot the tree we can see it is quite simple with 5 tips and 4 nodes:
    
    plot(tree, type = "cladogram", edge.width = 2)

Now, we will look at the structure of how R encodes phylogenies. This will help give a basic understanding of the right
path for phylogeny manipulation. First, we'll plot the tip and node labels on the tree:
    
    tiplabels()
    nodelabels()

Remember that the command "str" will look at the structure of an object in R. Let's do that with the phylogeny:
    str(tree)

We can see that our simple tree has three main components (edge, Nnode, and tip.label) and two attributes. We will look at
the three main components to see how they are structured:
    
    tree$Nnode
    tree$tip.label
    tree$edge

The Nnode and tip.label components should be pretty self explanatory. The edge table needs a little more explaining. In 
column 2, you have tips and nodes with their most immediate ancestor in column 1. So if you look at the plot of your simple
phylogeny with the tips and nodes labeled, you can see that the tips "A" and "B" are the tips labelled 1 and 2, and their 
most immediate ancestor is node 8. You can see this will match up with your edge table's values, where the rows with 1 and 2 
in the second column should have the number 8 in the first column. All in all, this edge table explains all of the 
relationships in the entire phylogenetic tree, and the larger the tree, the larger this table grows.

Now we'll do some modifications to our tree. You can do this in two ways, and we will do both. First, you can decide which 
tips you would like to keep in your pruned tree (and plot it to show the effect):
    
    tips_to_keep <- c("A", "B", "D", "E")
    
    pruned_tree1 <- keep.tip(tree, tips_to_keep)
    
    plot(pruned_tree1)
    
You can then check out the tip and node labels as well as the edge table. All of the labels of the tips and nodes and their
relationships will change if you add or remove tips from the phylogeny. 
  
    tiplabels()
    nodelabels()
    pruned_tree1$edge
    
You can also prune a tree by deciding which tips to remove. Here we will remove tip "D" and plot the tree again.
    
    tips_to_remove <- c("D")
    
    pruned_tree2 <- drop.tip(tree, tips_to_remove)
    
    plot(pruned_tree2)
    tiplabels()
    nodelabels()
    
    pruned_tree2$edge

Recall that we can rotate nodes on a tree, and the tree remains identical but visualized in a different way. Let's try 
rotating nodes on one of your pruned trees and plotting:
    
    # rotate at node 6 and node 7
    rotated_tree1 <- rotate(pruned_tree1, 6)
    
    plot(rotated_tree1)
    
    rotated_tree2 <- rotate(rotated_tree1, 7)
    
    plot(rotated_tree2)

One useful command to see if phylogenies are still identical after manipulation is the all.equal function. Let's try it here
with the double rotated tree with both of the pruned trees. Do you get the result you expect?
   
    all.equal(pruned_tree1, rotated_tree2)
    
    all.equal(pruned_tree2, rotated_tree2)
    
### 5. Let's work with a real phylogeny.
    
Let's read in the _Certhia_ tree and plot it.
    x <- read.nexus("certhia_summed.tre")
    plot(x, cex=0.7)

You can see that the tree is not rooted with the outgroup _C. familiaris_. Let's root it, we can do this two different ways.
    
    # root the tree with the outgroup
    x <- root(x, outgroup = "Cfamiliaris")
    plot(x, cex=0.7)

This rooted phylogeny doesn't look the prettiest with the outgroup outside the rest of the individuals. Let's try a 
different method:

    x <- midpoint.root(x)
    plot(x, cex=0.7)
    
Looks much better. Now let's prune this tree to just five individuals from Central America and plot:
    
    tips_to_keep <- c("Chiapas1", "Chiapas2", "Honduras1", "Honduras2", "Honduras3")
    
    x_trimmed <- keep.tip(x, tips_to_keep)
    
    plot(x_trimmed)

Now we'll measure the phylogenetic distance between tips. This will create a pairwise matrix of phylogenetic distances:

    cophenetic.phylo(x_trimmed)
    
The last thing we'll do, is write an output modified phylogenetic tree and read it back into R:

    write.tree(x_trimmed, file="test_tree.tre")
    
    test_tree <- read.tree("test_tree.tre")
    
    plot(test_tree)
