import random
import copy
import time
import os
import re
import operator
import deBruijnGraph

from utils import *
from test_utils import *
from DNAdroplet import DNADroplet
from DNAfountain import DNAFountain
from glass import Glass
# from deBruijnGraph import DeBruijnGraph


check_index = True
number_of_droplets = 10000
exp_seq_copy_num = 50
err_rate = 0.03
kmer_length = 18
kmer_cut = 4
data_block_length = 35
fountain_seed = 3
max_repeat_num_kmer = 5


primerF = 'CCTGCAGAGTAGCATGTC'  # 5'-->3'
primerE = 'CTGACACTGATGCATCCG'  # complement seq of P2

# file1 = open(r'../Apollo program.kux', 'rb')
work_dir = r'/work/0.Simulations/3.DecodingSpeedCompare/run' + str(fountain_seed) + '/'


# work_dir = r'E:\work\dec_speed' + '\\'
input_file = r'input_files/Dunhuang.6.8MB.zip'

# input_file = r'C:\Users\Lifu Song\Documents\CodingProjects\TinyCore.iso'

muscle = r'/home/lifu/Downloads/software/muscle3.8.31/muscle3.8.31_i86linux64 '  #Location of muscle program
DBGPS = r'/home/lifu/CodingProjects/switch-codes-switch-codes-kmer-cnt-DBG-master/DBGPS ' # Location of DBGPS-dy program
starcode = r'/home/lifu/starcode/starcode '


run_name = 'Speed Comparison'
file1 = open(input_file, 'rb')


fountain_init_index = 1
filebytes1 = file1.read()
file1.close()
# fdna1 = DNAFountain(filebytes1,data_block_length)
fdna1 = DNAFountain(filebytes1, data_block_length, fountain_init_index, fountain_seed)
fdna1.degree_table_folds = 2
fdna1.gen_degrees()

dps_seqs = {}
seq_arr = []

print('Generating Droplets')
droplets = get_droplets_check_repeat_kmer(number_of_droplets, fdna1, kmer_length)
for dps in droplets:
    dps_seqs[dps.head_index] = primerF + dps.to_DNA_CRC() + primerE
    seq_arr.append(primerF + dps.to_DNA_CRC() + primerE)
min_id = droplets[0].head_index
max_id = droplets[-1].head_index

dna_hd = DNAHandler()


print('Generating eDNAs')
eDNAs_all = []
for id in dps_seqs:
    eDNAs = dna_hd.copy_randnum([primerF + dps_seqs[id] + primerE ], exp_seq_copy_num, exp_seq_copy_num)
    eDNAs = dna_hd.add_rand_sub_new(eDNAs, err_rate / 2)
    eDNAs = dna_hd.add_rand_ins_new(eDNAs, err_rate / 4)
    eDNAs = dna_hd.add_rand_del_new(eDNAs, err_rate / 4)
    fo = open(work_dir + r'muscle/' + str(id), "tw")
    i = 1
    for seq in eDNAs:
        fo.write(">" + str(i) + "\n")
        fo.write(seq + "\n")
        i = i + 1
    fo.close()

    eDNAs_all.extend(eDNAs)

eDNAs_file_loc = work_dir + 'all_eDNAs.fa'
fo = open(eDNAs_file_loc,  "tw")

i = 1
for ee in eDNAs_all:
    fo.write(">" + str(i) + "\n")
    fo.write(ee + "\n")
    i = i + 1
fo.close()

res = {}
print('Recovering with DBGPS C version')
res['deG'] = {}
res['mus'] = {}

a = time.perf_counter()
print("Running DBGPS: ")
print(DBGPS + r' -k 21 -t 5 -c ' + str(kmer_cut) + r' -a ' + str(min_id) + r' -b ' + str(max_id) + r' ' + eDNAs_file_loc + r' >' + eDNAs_file_loc + r'.dec' + r' 2>' + eDNAs_file_loc + r'.dec.log' + str(fountain_seed))
os.system(DBGPS + r' -k 21 -t 5 -c ' + str(kmer_cut) + r' -a ' + str(min_id) + r' -b ' + str(max_id) + r' ' + eDNAs_file_loc + r' >' + eDNAs_file_loc + r'.dec' + r' 2>' + eDNAs_file_loc + r'.dec.log' + str(fountain_seed))

res['deG']['DBGPS'] = time.perf_counter() - a



print('Running starcode:')
a = time.perf_counter()
starcode_clu_file_loc = work_dir + r'starcode.clusters'
print(starcode + r'--print-clusters -t 5 -d 4 -s -i ' + eDNAs_file_loc + r' -o ' + starcode_clu_file_loc)
os.system(starcode + r'--print-clusters -t 5 -d 4 -s -i ' + eDNAs_file_loc + r' -o ' + starcode_clu_file_loc)
res['mus']['starcode'] = time.perf_counter() - a

print('Reading starcode clusters....')
clu_seqs = read_starcode_clusters(starcode_clu_file_loc, 3)
for clu in clu_seqs:
    clu_seqs_file = work_dir + 'muscle/' + str(clu)
    fo = open(clu_seqs_file, "tw")
    i=1
    for seq in clu_seqs[clu]:
        fo.write(">" + str(i) + "\n")
        fo.write(seq)
        fo.write("\n")
    fo.close()

a = time.perf_counter()
print('Running Muscle....')
for clu in clu_seqs:
    clu_seqs_file = work_dir + 'muscle/' + str(clu)
    os.system(muscle + r' -in ' + clu_seqs_file + r' -out ' + clu_seqs_file +  r'.aln')
res['mus']['muti-align'] = time.perf_counter() - a


a = time.perf_counter()
clu_cons_seqs = []
for clu in clu_seqs:
    clu_seqs_file = work_dir + 'muscle/' + str(clu)
    clu_aln_file = clu_seqs_file +  r'.aln'
    reads = read_fasta(clu_aln_file)
    cons = majority_merge(reads)
    clu_cons_seqs.append(cons)
    # if dps_seqs[id] in cons:
    #     i = i + 1
res['mus']['consensus'] = time.perf_counter() - a

dec_num = 0
for seq in dps_seqs:
    if seq_in_seqs(dps_seqs[seq], clu_cons_seqs):
        dec_num = dec_num + 1
res['mus']['Sr'] = dec_num





