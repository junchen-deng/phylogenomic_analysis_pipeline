# Analysis pipeline and post-processing scripts
Post-processing scripts for a single marker after running [**hybpiper retrieve_sequences**].  

**Mafft** or **fsa** were used to align aa sequences, which were reverse-translated to nucl alignment. **trimal** was used to trim gappy regions. The scripts also produce untrimmed alignment without outgroups, which is used to build trees with FastTree for the paralog filtering in the next step.  

![Phylogenomic_update_20240409](https://github.com/junchen-deng/post_processing_hybpiper/assets/75742791/6f0d8a4e-4394-4372-b617-2d7ef3c62820)

The R script **visualize_gene_trees.R** allowed a threshold-based method to narrow down the number of trees with potential long branches. Their DNA alignments were examined manually. The bar plot below shows the distribution of length ratios (Max. branch length/total branch length) from trees before the long-branch filtering step. We can see that most genes have ratios between 0 - 0.04. Then, we chose 0.03 as the threshold to be sure to include most genes in question. For the trees below the threshold, we took a quick look at their topology based on our knowledge of the relationships among sampled planthopper species. For other trees above the threshold, we manually examine their gene alignments.

![length_ratio_distribution](https://github.com/user-attachments/assets/4a3fcd71-b33d-4d61-949e-6bbc9a9c143b)
