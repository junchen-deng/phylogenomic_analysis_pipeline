#! /usr/bin/env python3

# requirement:
# mafft
# trimal
# FastTree

import os, glob, argparse, sys

parser = argparse.ArgumentParser(
			description='This script is a post-processing script for a single marker after running [hybpiper retrieve_sequences]. Mafft is applied to align aa sequences, which are revers-translated to nucl alignment. Then, trimal is used for trimming gappy regions. The script also produce untrimmed alignment without outgroups, which is used to build trees with FastTree for the paralog filtering in the next step.'	
)

# positional arguments
parser.add_argument("aa_path", metavar='<FASTA>', help="the path to the AA folder produced by [hybpiper retrieve_sequences]")
parser.add_argument("nucl_path", metavar='<FASTA>', help="the path to the nucl folder produced by [hybpiper retrieve_sequences]")
parser.add_argument("outdir", metavar='<FASTA>', type=str, help="output path")
parser.add_argument("gene_id", metavar='<str>', type=str, help="busco gene id")

args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

aa_path = args.aa_path
nucl_path = args.nucl_path
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
    

# 1) remove Xs and Ns in aa and nucl sequences <-- Xs and Ns were autoatically added by hybpiper and they will not be recognized by mafft  
# This step generates gene_id.FAA and gene_id.FAA.woOUT in the folder outdir/all_aa_woX
#                 and gene_id.FNA and gene_id.FNA.woOUT in the folder outdir/all_nucl_woN
id_outgroup1 = []
id_outgroup2 = []
for file in glob.glob("/home/junchen.deng/busco_test_fulgo/Homalodisca_vitripennis/*.faa"):
    label = file.split("/")[-1].split(".")[0]
    id_outgroup1.append(label)
for file in glob.glob("/home/junchen.deng/busco_test_fulgo/Callitettix_versicolor/*.faa"):
    label = file.split("/")[-1].split(".")[0]
    id_outgroup2.append(label)


prot_woX_woOUT = ImportFasta("{1}/{0}.FAA".format(gene_id, aa_path), 'aa')    # sequences without outgroups will be used to construct trees for filtering
prot_woX = list(prot_woX_woOUT)
nucl_woN_woOUT = ImportFasta("{1}/{0}.FNA".format(gene_id, nucl_path), 'nucl')
nucl_woN = list(nucl_woN_woOUT)
if gene_id in id_outgroup1:
    prot_outgroup1 = ImportFasta("/home/junchen.deng/busco_test_fulgo/Homalodisca_vitripennis/{}.faa".format(gene_id), 'aa')
    prot_outgroup1[0][0] = "Homalodisca_" + prot_outgroup1[0][0]    # give a specific header extension to each outgroup
    prot_woX += prot_outgroup1
    nucl_outgroup1 = ImportFasta("/home/junchen.deng/busco_test_fulgo/Homalodisca_vitripennis/{}.fna".format(gene_id), 'nucl')
    nucl_outgroup1[0][0] = "Homalodisca_" + nucl_outgroup1[0][0]
    nucl_woN += nucl_outgroup1
if gene_id in id_outgroup2:
    prot_outgroup2 = ImportFasta("/home/junchen.deng/busco_test_fulgo/Callitettix_versicolor/{}.faa".format(gene_id), 'aa')
    prot_outgroup2[0][0] = "Callitettix_" + prot_outgroup2[0][0]
    prot_woX += prot_outgroup2
    nucl_outgroup2 = ImportFasta("/home/junchen.deng/busco_test_fulgo/Callitettix_versicolor/{}.fna".format(gene_id), 'nucl')
    nucl_outgroup2[0][0] = "Callitettix_" + nucl_outgroup2[0][0]
    nucl_woN += nucl_outgroup2        
with open("{1}/all_aa_woX/{0}.FAA".format(gene_id, outdir), "w") as file_out:
    for entry in prot_woX: 
        print('>', entry[0], '\n', entry[1], '\n', end = '', sep = '', file = file_out)
with open("{1}/all_aa_woX/{0}.FAA.woOUT".format(gene_id, outdir), "w") as file_out:
    for entry in prot_woX_woOUT: 
        print('>', entry[0], '\n', entry[1], '\n', end = '', sep = '', file = file_out)
with open("{1}/all_nucl_woN/{0}.FNA".format(gene_id, outdir), "w") as file_out:
    for entry in nucl_woN: 
        print('>', entry[0], '\n', entry[1], '\n', end = '', sep = '', file = file_out)
with open("{1}/all_nucl_woN/{0}.FNA.woOUT".format(gene_id, outdir), "w") as file_out:
    for entry in nucl_woN_woOUT: 
        print('>', entry[0], '\n', entry[1], '\n', end = '', sep = '', file = file_out)         

# 2) generate aa alignment with mafft and reverse translate aa alignment to nucl alignment
# stop_codon = ['TAA', 'TAG', 'TGA']; mafft delete stop codon "*" in aa alginment but there are stop codons in nucl seqs 
# This step generates gene_id.align and gene_id.align.woOUT in the folder outdir/all_aa_align
#                 and gene_id.align and gene_id.align.woOUT in the folder outdir/all_nucl_align

os.system("mafft --thread 6 --maxiterate 1000 --quiet --localpair {1}/all_aa_woX/{0}.FAA > {1}/all_aa_align/{0}.align".format(gene_id, outdir))
os.system("mafft --thread 6 --maxiterate 1000 --quiet --localpair {1}/all_aa_woX/{0}.FAA.woOUT > {1}/all_aa_align/{0}.align.woOUT".format(gene_id, outdir))

nucl_woN = ImportFasta("{1}/all_nucl_woN/{0}.FNA".format(gene_id, outdir), 'nucl')
nucl_woN_woOUT = ImportFasta("{1}/all_nucl_woN/{0}.FNA.woOUT".format(gene_id, outdir), 'nucl')
aa_align = ImportFasta("{1}/all_aa_align/{0}.align".format(gene_id, outdir), 'aa')
aa_align_woOUT = ImportFasta("{1}/all_aa_align/{0}.align.woOUT".format(gene_id, outdir), 'aa')

nucl_align = RevTrans(aa_align, nucl_woN)
nucl_align_woOUT = RevTrans(aa_align_woOUT, nucl_woN_woOUT)

with open("{1}/all_nucl_align/{0}.align".format(gene_id, outdir), "w") as file_out:
    for entry in nucl_align: 
        print('>', entry[0], '\n', entry[1], '\n', end = '', sep = '', file = file_out)
with open("{1}/all_nucl_align/{0}.align.woOUT".format(gene_id, outdir), "w") as file_out:
    for entry in nucl_align_woOUT: 
        print('>', entry[0], '\n', entry[1], '\n', end = '', sep = '', file = file_out)  

# 3) generate trees from nucl alignment without out groups using FastTree
# This step generates gene_id.paralogs.tre in outdir/all_treefile
os.system("FastTreeMP -nt -gtr {1}/all_nucl_align/{0}.align.woOUT > {1}/all_treefile/{0}.paralogs.tre".format(gene_id, outdir))

# 4) trim aa alignment with outgroups using trimal; then produce trimmed nucl alignment
# !!! trimming will remove any zero-length seqs; during reverse translation, the correpsonding sample will be empty with only a title
os.system("trimal -in {1}/all_aa_align/{0}.align -out {1}/all_aa_align/{0}.trim -automated1 -colnumbering > {1}/all_aa_align/{0}.txt".format(gene_id, outdir)) 

nucl_align = ImportFasta("{1}/all_nucl_align/{0}.align".format(gene_id, outdir), 'nucl')
pos_file = open("{1}/all_aa_align/{0}.txt".format(gene_id, outdir), 'r')
pos = pos_file.read()
pos_list = pos.split(', ')
pos_list[-1] = pos_list[-1].replace('\n','')
pos_list[0] = pos_list[0].replace('#ColumnsMap\t','')    # pos_list = ['3','4','7'...]

nucl_trim = []
for sample in nucl_align:
    trim = ''
    for item in pos_list: 
        i = int(item)
        keep = sample[1][i*3 : (i+1)*3]
        trim += keep
    nucl_trim.append([sample[0],trim])

with open("{1}/all_nucl_align/{0}.trim".format(gene_id, outdir), "w") as file_out:          
    for entry in nucl_trim:
        print('>', entry[0], '\n', entry[1], '\n', end = '', sep = '', file = file_out)      


        
     




