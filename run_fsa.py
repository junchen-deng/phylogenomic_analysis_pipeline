#! /usr/bin/env python3

# requirement:
# mafft
# trimal
# FastTree

import os, glob, argparse, sys

parser = argparse.ArgumentParser(
			description='This script is a post-processing script for a single marker after running [hybpiper retrieve_sequences]. fsa is applied to align aa sequences, which are revers-translated to nucl alignment. Then, trimal is used for trimming gappy regions. The script also produce untrimmed alignment without outgroups, which is used to build trees with FastTree for the paralog filtering in the next step.'	
)

# positional arguments
# parser.add_argument("aa_path", metavar='<FASTA>', help="the path to the AA folder produced by [hybpiper retrieve_sequences]")
# parser.add_argument("nucl_path", metavar='<FASTA>', help="the path to the nucl folder produced by [hybpiper retrieve_sequences]")
parser.add_argument("outdir", metavar='<FASTA>', type=str, help="output path")
parser.add_argument("gene_id", metavar='<str>', type=str, help="busco gene id")

args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

# aa_path = args.aa_path
# nucl_path = args.nucl_path
outdir = args.outdir
gene_id = args.gene_id

def ImportFasta(fasta_file, data_type):
   FASTA = open(fasta_file, 'r')
   Seq_list = []
   Sequence = ''
   Seq_heading = ''
   for line in FASTA:   # Copying the sequence (potentially spread across multiple lines) to a single line
      if line.startswith('>'):
         if Sequence != '':   # Saves the existing Seq_heading and Sequence to a list before overwriting them
            if data_type == 'aa':
                Sequence = Sequence.replace('X','')
            elif data_type == 'nucl':
                Sequence = Sequence.replace('N','')
            Seq_list.append([Seq_heading, Sequence])
         Sequence = ''
         Seq_heading = line.strip().strip(">") # Takes the whole name as heading
      else:
         Sequence = Sequence + line.strip().upper()
   if data_type == 'aa':
      Sequence = Sequence.replace('X','')
   elif data_type == 'nucl':
      Sequence = Sequence.replace('N','')
   Seq_list.append([Seq_heading, Sequence.strip('')]) # Saves the final sequence (Seq_heading and Sequence) to a list
   FASTA.close()
   return(Seq_list)

def RevTrans(aa_align, nucl):    # should be both list of list
    stop_codon = ['TAA', 'TAG', 'TGA']
    nucl_align = []
    for i in range(len(aa_align)):
        header = aa_align[i][0]
        nucl_align.append([header,''])
        aa_seq = aa_align[i][1]
        nucl_seq = nucl[i][1]
        
        nucl_seq_clean = ''    # remove potential stop codons captured by hybpiper
        x = 1
        for k in range(int(len(nucl_seq)/3)):
            if nucl_seq[(x-1)*3 : x*3] not in stop_codon:   
                nucl_seq_clean += nucl_seq[(x-1)*3 : x*3]
            x += 1
            
        nucl_seq = nucl_seq_clean    # now nucl_seq is the same length as aa_seq
        x = 1
        for j in range(len(aa_seq)):
            if aa_seq[j]  != '-':
                nucl_align[i][1] += nucl_seq[(x-1)*3 : x*3]                    
            else:
                nucl_align[i][1] += '---'
                x -= 1
            x += 1        
    return(nucl_align)
    
# fsa alignment
os.system("fsa --fast {1}/all_aa_align/{0}.trim > {1}/all_aa_fsa/{0}.align".format(gene_id, outdir))

# remove outgroups in xxx.align; prepare .woOUT for FastTree
remove_list = ['Callitettix','Homalodisca']
prot_list = ImportFasta("{1}/all_aa_fsa/{0}.align".format(gene_id, outdir),'aa')
prot_list_new = list(prot_list)
for gene in prot_list:
    if gene[0].split('_')[0] in remove_list:
        prot_list_new.remove(gene)
if prot_list_new == []:
    print("{0} is empty after removing the marker".format(gene_id))
with open("{1}/all_aa_fsa/{0}.align.woOUT".format(gene_id, outdir),'w') as file_out: 
    for entry in prot_list_new:
        print('>', entry[0], '\n', entry[1], '\n', end = '', sep = '', file = file_out)    

nucl_woN = ImportFasta("{1}/all_nucl_align/{0}.trim".format(gene_id, outdir), 'nucl')
nucl_woN_new = list(nucl_woN)
for gene in nucl_woN:    # create unaligned seqs
    gene[1] = gene[1].replace('-','')

nucl_woN_new = list(nucl_woN)    # removing empty sequences from unaligned sequence!!!
for gene in nucl_woN:
    if gene[1] == '':
        nucl_woN_new.remove(gene)    
        
aa_align = ImportFasta("{1}/all_aa_fsa/{0}.align".format(gene_id, outdir), 'aa')
nucl_align = RevTrans(aa_align, nucl_woN_new)
nucl_align_woOUT = list(nucl_align)
# remove outgroups from nucl
for gene in nucl_align:
    if gene[0].split('_')[0] in remove_list:
        nucl_align_woOUT.remove(gene)
        
with open("{1}/all_nucl_fsa/{0}.align".format(gene_id, outdir), "w") as file_out:
    for entry in nucl_align: 
        print('>', entry[0], '\n', entry[1], '\n', end = '', sep = '', file = file_out)
with open("{1}/all_nucl_fsa/{0}.align.woOUT".format(gene_id, outdir),'w') as file_out: 
    for entry in nucl_align_woOUT:
        print('>', entry[0], '\n', entry[1], '\n', end = '', sep = '', file = file_out)  

# generate trees from nucl alignment without out groups using FastTree
# This step generates gene_id.paralogs.tre in outdir/all_treefile
os.system("FastTreeMP -nt -gtr {1}/all_nucl_fsa/{0}.align.woOUT > {1}/all_treefile/{0}.paralogs.tre".format(gene_id, outdir))

# trim aa alignment with outgroups using trimal; then produce trimmed nucl alignment
os.system("trimal -in {1}/all_aa_fsa/{0}.align -out {1}/all_aa_fsa/{0}.trim -automated1 -colnumbering > {1}/all_aa_fsa/{0}.txt".format(gene_id, outdir)) 

pos_file = open("{1}/all_aa_fsa/{0}.txt".format(gene_id, outdir), 'r')
pos = pos_file.read()
pos_list = pos.split(', ')
pos_list[-1] = pos_list[-1].replace('\n','')
pos_list[0] = pos_list[0].replace('#ColumnsMap\t','')    # pos_list = ['3','4','7'...]

nucl_trim = []
for sample in nucl_align:    # trimal will remove any empty sequences in aa_trim; "nucl_align" can contain empty sequences by trimming in this way!!!! 
    trim = ''
    for item in pos_list: 
        i = int(item)
        keep = sample[1][i*3 : (i+1)*3]
        trim += keep
    nucl_trim.append([sample[0],trim])

with open("{1}/all_nucl_fsa/{0}.trim".format(gene_id, outdir), "w") as file_out:          
    for entry in nucl_trim:
        print('>', entry[0], '\n', entry[1], '\n', end = '', sep = '', file = file_out)      


        
     




