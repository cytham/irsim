"""

Execute.py

This module simulates FASTQ short-reads using DWGSIM.

Copyright (C) 2019 Tham Cheng Yong

This file is part of IRSim.

IRSim is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

IRSim is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with IRSim.  If not, see <https://www.gnu.org/licenses/>.
"""

import os.path
import numpy as np
import multiprocessing
from collections import OrderedDict
from subprocess import DEVNULL, STDOUT, check_call

#Decipher amount of simulated reads to suit x% of transcripts FPKM
def sim_amount(gene_list, fpkm_dict, num_read, total_len, percentile_val, fasta_dict, total_reps):
    #Get xth Percentile FPKM
    FPKM = []
    for gene in gene_list:
        FPKM.append(float(fpkm_dict[gene]))
    percentile_FPKM = np.percentile(np.array(FPKM),percentile_val)
    #Average fragment depth = (number of fragments x fragment size)/Total bases
    #Number of fragments per gene = average fragment depth x (cDNA length/fragment size)
    #FPKM = (Number of fragments per gene/(Number of reads/1000000))/(cDNA length/1000)
    #Merge and simplify
    #FPKM = (Number of fragments x 1000000000)/(Total bases x number of reads)
    #Number of fragments = (FPKM x Total bases x total num of fragments)/1000000000
    #Adjust "Total bases" by ratio of total_transcripts/number of genes
    #Number of fragments = (FPKM x Total bases x total num of fragments)/(1000000000) multiply by number of replicates  # x (total_transcripts/number of genes)
    sim1_num_frags = int((percentile_FPKM*total_len*num_read/2)/(1000000000))*total_reps  #*(len(fasta_dict)/len(gene_list))
    return sim1_num_frags, percentile_FPKM

#Decipher amount of simulated reads to suit x% of transcripts fragments
def sim_amount2(gene_list, fpkm_dict, total_len, percentile_val, total_reps, transcript_frag_dict, transcript_dict, cdna_size_dict):
    #Get xth Percentile FPKM
    FPKM_dict = OrderedDict()
    FPKM = []
    for gene in gene_list:
        FPKM_dict[gene] = float(fpkm_dict[gene])
    FPKM_tup_sort = sorted(FPKM_dict.items(), key=lambda x:x[1])
    for gene,fpkm in FPKM_tup_sort:
        FPKM.append(fpkm)
    percentile_index = int(len(FPKM)*(percentile_val/100)) - 1
    gene_name, percentile_FPKM = FPKM_tup_sort[percentile_index]
    for rep in transcript_frag_dict:
        last_rep_name = rep
    #Average fragment depth = (number of fragments x fragment size)/Total bases
    #Number of fragments per gene = average fragment depth x (cDNA length/fragment size)
    #FPKM = (Number of fragments per gene/(Number of reads/1000000))/(cDNA length/1000)
    #Merge and simplify
    #FPKM = (Number of fragments x 1000000000)/(Total bases x number of reads)
    #Number of fragments = (FPKM x Total bases x number of reads)/1000000000
    #Adjust "Total bases" by ratio of total_transcripts/number of genes
    #Number of fragments = (FPKM x Total bases x number of reads)/(1000000000)) multiply by number of replicates
    #Number of fragments = (Number of fragments per gene x Total bases)/cDNA length) multiply by number of replicates
    try:
        sim1_num_frags = int((transcript_frag_dict[last_rep_name]['@'+transcript_dict[gene_name]+'-e']*total_len)/cdna_size_dict[transcript_dict[gene_name]])*total_reps
    except:
        sim1_num_frags = int((transcript_frag_dict[last_rep_name]['@'+transcript_dict[gene_name]+'-i']*total_len)/cdna_size_dict[transcript_dict[gene_name]])*total_reps
    return sim1_num_frags, percentile_FPKM

#Chunk list generator
def split_chunks(a, n):
    k, m = divmod(len(a), n)
    return list(a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))

#Split reference fasta to multiple fasta files according to number of threads
def split_fasta(output_dir, ref, num_threads, num_fasta):
    fasta = os.path.join(output_dir, ref)
    tmp_dict = OrderedDict()
    tmp_list = []
    for index in range(num_threads):
        tmp_dict[index] = open(os.path.join(output_dir, 'tmp' + str(index) + '.fa'), 'w')
    with open(fasta) as f:
        for line in f:
            tmp_list.append(line + next(f))
    #Split into chunks
    chunk_list = split_chunks(range(num_fasta), num_threads)
    for index in range(len(chunk_list)):
        for num in chunk_list[index]:
            tmp_dict[index].write(tmp_list[num])
    for index in tmp_dict:
        tmp_dict[index].close()

#Run single dwgsim
def single_subprocess_dwgsim(dwgsim_path, output_dir, num_threads, first_read_error, second_read_error, ins_length, ins_length_stdev, split_frags, read_length, mut_rate, indel_rate, indel_ext, indel_len, ran_read, seed, ref):
    check_call([str(dwgsim_path), '-e', str(first_read_error), '-E', str(second_read_error), '-d', str(ins_length), '-s', str(ins_length_stdev), '-N', str(split_frags), '-1', str(read_length), '-2', str(read_length), '-r', str(mut_rate), '-R', str(indel_rate), '-X', str(indel_ext), '-I', str(indel_len), '-y', str(ran_read), '-z', str(seed), '-o', '1', '-c', '0', '-S', '2', str(ref), str(ref + '.out')], stdout=DEVNULL, stderr=STDOUT)

#Run multi dwgsim
def multi_subprocess_dwgsim(current_gene_list, dwgsim_dir, output_dir, num_threads, first_read_error, second_read_error, ins_length, ins_length_stdev, sim_num_frags, read_length, mut_rate, indel_rate, indel_ext, indel_len, ran_read, seed):
    dwgsim_path = os.path.join(dwgsim_dir, 'dwgsim')
    in_fasta_list = []
    if len(current_gene_list) >= num_threads:
        for index in range(num_threads):
            in_fasta_list.append(os.path.join(output_dir, 'tmp' + str(index) + '.fa'))
        split_frags = round(sim_num_frags/num_threads)
        jobs = []
        for ref in in_fasta_list:
            p = multiprocessing.Process(target=single_subprocess_dwgsim, args=(dwgsim_path, output_dir, num_threads, first_read_error, second_read_error, ins_length, ins_length_stdev, split_frags, read_length, mut_rate, indel_rate, indel_ext, indel_len, ran_read, seed, ref))
            jobs.append(p)
            p.start()
    else:
        for index in range(len(current_gene_list)):
            in_fasta_list.append(os.path.join(output_dir, 'tmp' + str(index) + '.fa'))
        split_frags = round(sim_num_frags/len(current_gene_list))
        jobs = []
        for ref in in_fasta_list:
            p = multiprocessing.Process(target=single_subprocess_dwgsim, args=(dwgsim_path, output_dir, num_threads, first_read_error, second_read_error, ins_length, ins_length_stdev, split_frags, read_length, mut_rate, indel_rate, indel_ext, indel_len, ran_read, seed, ref))
            jobs.append(p)
            p.start()
    for job in jobs:
        job.join()
    #for ref in in_fasta_list:
        #check_call(['rm', ref], stdout=DEVNULL, stderr=STDOUT)
        #check_call(['rm', ref + '.out.mutations.vcf'], stdout=DEVNULL, stderr=STDOUT)
        #check_call(['rm', ref + '.out.mutations.txt'], stdout=DEVNULL, stderr=STDOUT)
