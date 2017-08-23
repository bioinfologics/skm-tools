# Working with skip-mers

This repository contains a number of tools that work with skip-mers.

# The concept of skip-mer

Skip-mers are sequences of nucleotides taken in cycles of m every n for a total of k. The full description is in the biorXiv manuscript: [Skip-mers: increasing entropy and sensitivity to detect conserved genic regions with simple cyclic q-grams](http://www.biorxiv.org/content/early/2017/08/23/179960)

It is actually a lot easier to show with a picture:

![Definition image](definition.png?raw=true)

From the manuscript:

>A skip-mer is a simple cyclic q-gram that includes m out of every n bases until a total of k bases is reached. Its shape is defined by a function SkipMer(m,n,k), as shown in Figure 1. To maintain cyclic properties and the existence of the reverse-complement as a skip-mer defined by the same function, k must be a multiple of m. This also enables the existence of a canonical representation for each skip-mer, defined as the lexicographically smaller of the forward and reverse-complement representations.


# Tools used in the manuscript

## skm-count

This is a basic skip-mer counter inputs are sequences, output is a skipmer spectrum file (binary) and a text skip-mer
spectrum histogram, with a maximum unique skip-mer progression plot.


## skm-multiway-coverage

This creates a multiway coverage analysis of a reference sequence vs. a set of alternative sequences. In the publication
it is used to showcase how gene analysis through skip-mers is more robust than with kmers when
including multiple samples across a phylogeny.

The tool loads a reference genome and a set of features from a gff3 file, alongside a set of "coverage genomes". Then it creates a SkipMerPosition index, and projects the skip-mer coverage of all the other datasets into it. The coverage form the progressive incorporation of datasets is dumped, and a sequence conservation score for each feature versus each  "coverage genome". The "detailed" scores are based on skip-mer counts, and the "projected" are based on covered bases as described in the manuscript.


# Tools in development (use them at your own risk)

More tools are comming soon, you can have a look in the alternative branches, but no warranty given.
