# Analysis pipeline and post-processing scripts
Post-processing scripts for a single marker after running [**hybpiper retrieve_sequences**].  
**Mafft** or **fsa** were used to align aa sequences, which were reverse-translated to nucl alignment. **trimal** was used to trim gappy regions.  
The scripts also produce untrimmed alignment without outgroups, which is used to build trees with FastTree for the paralog filtering in the next step.  

The R script **visualize_gene_trees.R** allowed a threshold-based method to narrow down the number of trees with potential long branches. Their DNA alignments were examined manually.  
**length_ratio_distribution.png** showed the distribution of length ratios (Max. branch length/total branch length) from trees before the long-branch filtering step. 


![Phylogenomic_update_20240409](https://github.com/junchen-deng/post_processing_hybpiper/assets/75742791/6f0d8a4e-4394-4372-b617-2d7ef3c62820)
