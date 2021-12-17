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


work_dir = r'/work/1.DeBruijnGraphDecoding/6M/1.rawdata/'


sim_file = work_dir + r'6.5MB.DNAs.newids.tab.nohead'
input_file = work_dir + r'P10_5_BDDP210000009-1A_join_pear.fq.assembled.fastq.cluster.fq.t10.d4'

muscle = r'/home/lifu/software/muscle3.8.31/muscle3.8.31_i86linux64 '  #Location of muscle program

os.system('mkdir ' + work_dir + 'muscle')

max_seq_copy = 100000

run_name = 'Speed Comparison'



dps_seqs = read_sim(sim_file)

res = {}
res['muti-align'] = []
print('Reading starcode clusters....')
clu_seqs = read_starcode_clusters(input_file, 3)
for clu in clu_seqs:
    clu_seqs_file = work_dir + 'muscle/' + str(clu)
    fo = open(clu_seqs_file, "tw")
    i=1
    for seq in clu_seqs[clu]:
        if i > max_seq_copy:
            break
        fo.write(">" + str(i) + "\n")
        fo.write(seq)
        fo.write("\n")
        i = i + 1
    fo.close()

a = time.perf_counter()
print('Running Muscle....')
for clu in clu_seqs:
    clu_seqs_file = work_dir + 'muscle/' + str(clu)
    os.system(muscle + r' -align ' + clu_seqs_file + r' -output ' + clu_seqs_file + r'.aln')
    #os.system(muscle + r' -in ' + clu_seqs_file + r' -out ' + clu_seqs_file +  r'.aln')
    #print(muscle + r' -in ' + clu_seqs_file + r' -out ' + clu_seqs_file +  r'.aln')
res['muti-align'].append(time.perf_counter() - a)

a = time.perf_counter()
print('Running Muscle....')
for clu in clu_seqs:
    clu_seqs_file = work_dir + 'muscle/' + str(clu)
    #os.system(muscle + r' -align ' + clu_seqs_file + r' -output ' + clu_seqs_file + r'.aln') #muscle5
    os.system(muscle + r' -in ' + clu_seqs_file + r' -out ' + clu_seqs_file +  r'.aln') # muscle 3
res['muti-align'].append(time.perf_counter() - a)

a = time.perf_counter()
print('Running Muscle....')
for clu in clu_seqs:
    clu_seqs_file = work_dir + 'muscle/' + str(clu)
    os.system(muscle + r' -align ' + clu_seqs_file + r' -output ' + clu_seqs_file + r'.aln')
    #os.system(muscle + r' -in ' + clu_seqs_file + r' -out ' + clu_seqs_file +  r'.aln')
res['muti-align'].append(time.perf_counter() - a)


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
res['consensus'].append(time.perf_counter() - a)


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
res['consensus'].append(time.perf_counter() - a)


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
res['consensus'].append(time.perf_counter() - a)

dec_num = check_cons(clu_cons_seqs, dps_seqs)
#for seq in dps_seqs:
#    if seq_in_seqs(dps_seqs[seq], clu_cons_seqs):
#        dec_num = dec_num + 1
res['Sr'] = dec_num



a = time.perf_counter()
print('Running Muscle with perfect clusters....')
for id in dps_seqs:
    id_seqs_file = work_dir + 'muscle/' + str(id)
    os.system(muscle + r' -in ' + id_seqs_file + r' -out ' + id_seqs_file +  r'.aln')
res['mus']['muti-align-perfect-clusters'] = time.perf_counter() - a


a = time.perf_counter()
clu_cons_seqs = []
for id in dps_seqs:
    clu_seqs_file = work_dir + 'muscle/' + str(id)
    clu_aln_file = clu_seqs_file +  r'.aln'
    reads = read_fasta(clu_aln_file)
    cons = majority_merge(reads)
    clu_cons_seqs.append(cons)
    # if dps_seqs[id] in cons:
    #     i = i + 1
res['mus']['consensus-perfect-clusters'] = time.perf_counter() - a

dec_num = 0
for seq in dps_seqs:
    if seq_in_seqs(dps_seqs[seq], clu_cons_seqs):
        dec_num = dec_num + 1
res['mus']['Sr-perfect-clusters'] = dec_num



# io_time = []
#
# a = time.perf_counter()
# for id in range(0, 10000):
#     id_seq_file = work_dir + '1seq'
#     os.system(muscle + r' -in ' + id_seq_file + r' -out ' + id_seq_file +  r'.aln')
#
# io_time.append(time.perf_counter() - a)
#
#
# a = time.perf_counter()
# for id in range(0, 10000):
#     id_seq_file = work_dir + '1seq'
#     os.system(muscle + r' -in ' + id_seq_file + r' -out ' + id_seq_file +  r'.aln')
#
# io_time.append(time.perf_counter() - a)
#
#
# a = time.perf_counter()
# for id in range(0, 10000):
#     id_seq_file = work_dir + '1seq'
#     os.system(muscle + r' -in ' + id_seq_file + r' -out ' + id_seq_file +  r'.aln')
#
# io_time.append(time.perf_counter() - a)


