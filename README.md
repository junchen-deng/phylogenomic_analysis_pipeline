# post_processing_hybpiper
post-processing scripts for a single marker after running [hybpiper retrieve_sequences]. 
Mafft or fsa is applied to align aa sequences, which are reverse-translated to nucl alignment. Then, trimal is used for trimming gappy regions. 
The scripts also produce untrimmed alignment without outgroups, which is used to build trees with FastTree for the paralog filtering in the next step.
