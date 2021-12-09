from test_utils import *

kmer_length = 27
data_block_length = 35
fountain_seed = 2

p1 = 'CCTGCAGAGTAGCATGTC'  # 5'-->3'
p2 = 'CTGACACTGATGCATCCG'  # complement seq of P2

work_dir = r'input_files/'
file1 = open(work_dir + r'Dunhuang.6.8MB.zip', 'rb')

fountain_init_index = 101010101
filebytes1 = file1.read()
file1.close()
fdna1 = DNAFountain(filebytes1, data_block_length, fountain_init_index, fountain_seed, 16, 16, 8)
fdna1.ROBUST_FAILURE_PROBABILITY = 0.01
fdna1.c_value = 0.01
fdna1.gen_degrees()

dna_hd = DNAHandler()
i = 0
res = {}

for strand_num in [210000]:
    res[strand_num] = {}
    droplet_all = get_droplets_check_repeat_kmer(strand_num, fdna1, kmer_length)

    dnas = []
    for dp in droplet_all:
        dnas.append(p1 + dp.to_DNA_CRC_sIndex() + p2)

    for cov in [10]:
        res[strand_num][cov] = {}
        err_rate = 0.001
        err_rate_step = 0.001

        while err_rate <= 0.1:
            print(err_rate)
            eDNAs = dna_hd.copy_seqs(dnas, cov)
            eDNAs = dna_hd.add_rand_sub_new(eDNAs, err_rate/2, True)
            eDNAs = dna_hd.add_rand_indel_new(eDNAs, err_rate/4, err_rate/4, True)
            res[strand_num][cov][err_rate] = kdkn(eDNAs, dnas, kmer_length)
            err_rate = err_rate + err_rate_step

