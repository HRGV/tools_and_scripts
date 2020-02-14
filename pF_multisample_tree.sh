#!/bin/bash
#script to collect the SPADES assembled rRNAS and references in a common tree

#variables to be set
project=SYLT19MM

#Extract all PF spades assembled 16S rRNA gene sequencens and the related SILVA dbhits (and possibly also emirge hits)
for i in *tar.gz; do  tar -xf $i ${i%phy*}spades_rRNAs.final.fasta; done
for i in *tar.gz; do  tar -xf $i ${i%phy*}all.dbhits.NR97.fa; done

#concatenate PF hits and dbhits
cat *.spades_rRNAs.final.fasta > $project.collection_spades.fasta
cat *.all.dbhits.NR97.fa > $project.collection_references.fasta

#run the dbits through vsearch to remove duplicates 
vsearch --cluster_fast $project.collection_references.fasta --id 0.97 --relabel_keep --notrunclabels --centroids $project.collection_references_NR97.fasta

#concatenate PF hits and references
cat $project.collection_spades.fasta $project.collection_references_NR97.fasta > $project.spades_and_refsNR97_collection.fasta

#calculate alignment 
mafft-xinsi --thread 24 --adjustdirection $project.spades_and_refsNR97_collection.fasta > $project.spades_and_refsNR97_collection.mafft_xinsi.fasta

#calculate tree
fasttree -gtr -nt -quote $project.spades_and_refsNR97_collection.mafft_xinsi.fasta > $project.spades_and_refsNR97_collection.mafft_xinsi.fasttree_GTR.tre
