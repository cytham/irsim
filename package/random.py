"""
Random.py

This module creates random reads.

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

import gzip
import pickle
import os.path
import multiprocessing
from package.read import take_read_id
from collections import OrderedDict

#Count number of reads
def read_count(output_dir, rep_names_dict):
    read_count_dict = OrderedDict()
    for rep in rep_names_dict:
        read1_path = os.path.join(output_dir, 'sim_' + str(rep_names_dict[rep]) + '_read1.fastq')
        line_counts = 0
        with open(read1_path) as f:
            for line in f:
                line_counts += 1
        read_count_dict[rep] = line_counts/4
    return dict(read_count_dict)

#Random read generation
def random_read_counter(read_count_dict, num_read, ran_read):
    num_ran_reads = 0
    for rep in read_count_dict:
        if read_count_dict[rep] >= num_read*(1-ran_read)/2:
            read_count_dict[rep] = num_read*(1-ran_read)/2
        else:
            pass
        num_ran_reads = num_ran_reads + (num_read/2 - read_count_dict[rep])
    return int(num_ran_reads) 

#Create mock fasta for random read generation
def mock_fasta(output_dir, num_threads):
    name = '>mock_sequence'
    motif = 'GCAT'
    sequence = motif*1000
    for i in range(num_threads):
        file_path = os.path.join(output_dir, 'tmp' + str(i) + '.fa')
        out = open(file_path, 'w')
        out.write('{}\n{}\n'.format(name, sequence))
        out.close()

#Extract random fastq (Multi_thread, Pickle)
def extract_random_fastq_pickle(fastq_name, output_dir, trans_gene_dict):
    fastq1 = os.path.join(output_dir, fastq_name + '.read1.fastq.gz')
    fastq2 = os.path.join(output_dir, fastq_name + '.read2.fastq.gz')
    fastq_dict1 = OrderedDict()
    fastq_dict2 = OrderedDict()
    read_name_list = []
    with gzip.open(fastq1, 'rt') as f1:
        for line in f1:
            if line.rsplit('_', 9)[0] not in read_name_list:
                read_name_list.append(line.rsplit('_', 9)[0])
            try:
                fastq_dict1[line.rsplit('_', 9)[0]].append([line, next(f1), next(f1), next(f1)])
            except:
                fastq_dict1[line.rsplit('_', 9)[0]] = []
                fastq_dict1[line.rsplit('_', 9)[0]].append([line, next(f1), next(f1), next(f1)])
    with gzip.open(fastq2, 'rt') as f2:
        for line in f2:
            try:
                fastq_dict2[line.rsplit('_', 9)[0]].append([line, next(f2), next(f2), next(f2)])
            except:
                fastq_dict2[line.rsplit('_', 9)[0]] = []
                fastq_dict2[line.rsplit('_', 9)[0]].append([line, next(f2), next(f2), next(f2)])
    for name in read_name_list:
        fastq_dict1[name].sort(key=take_read_id)
        fastq_dict2[name].sort(key=take_read_id)
    pickle1 = open(os.path.join(output_dir, fastq_name) + '1.pkl', 'wb')
    pickle2 = open(os.path.join(output_dir, fastq_name) + '2.pkl', 'wb')
    pickle.dump(fastq_dict1, pickle1, -1)
    pickle.dump(fastq_dict2, pickle2, -1)
    pickle1.close()
    pickle2.close()

#Multiprocessing extract random fastq (Pickle)
def multi_thread_extract_random_fastq_pickle(output_dir, trans_gene_dict, strand_dict, num_threads):    
    #manager = multiprocessing.Manager()
    #dict1 = manager.dict()
    #dict2 = manager.dict()
    jobs = []
    for index in range(int(num_threads)):
        p = multiprocessing.Process(target=extract_random_fastq_pickle, args=('tmp' + str(index) + '.fa.out.bwa', output_dir, trans_gene_dict))
        jobs.append(p)
        p.start()
    for job in jobs:
        job.join()
    fastq_dict1 = OrderedDict()
    fastq_dict2 = OrderedDict()
    fastq_dict1['@rand'] = []
    fastq_dict2['@rand'] = []
    #Unpickle
    count_1 = 0
    count_2 = 0
    for index in range(int(num_threads)):
        pickle1 = open(os.path.join(output_dir, 'tmp' + str(index) + '.fa.out.bwa') + '1.pkl', 'rb')
        pickle2 = open(os.path.join(output_dir, 'tmp' + str(index) + '.fa.out.bwa') + '2.pkl', 'rb')
        dict1 = pickle.load(pickle1)
        dict2 = pickle.load(pickle2)
        pickle1.close()
        pickle2.close()
        for read in dict1['@rand']:
            fastq_dict1['@rand'].append(['@rand_0_0_0_0_1_1_0:0:0_0:0:0_' + str(count_1) + '/1\n'] + read[1:])
            count_1 += 1
        for read in dict2['@rand']:
            fastq_dict2['@rand'].append(['@rand_0_0_0_0_1_1_0:0:0_0:0:0_' + str(count_2) + '/2\n'] + read[1:])
            count_2 += 1
    return dict(fastq_dict1), dict(fastq_dict2)

#Write random reads
def random_reads_write(read_count_dict, ran_fastq_dict1, ran_fastq_dict2, output_dir, rep_names_dict, num_read):
    counter = 0
    for rep in rep_names_dict:
        fastq_list1 = ran_fastq_dict1['@rand'][counter:int(counter+(num_read/2 - int(read_count_dict[rep])))]
        fastq_list2 = ran_fastq_dict2['@rand'][counter:int(counter+(num_read/2 - int(read_count_dict[rep])))]
        write_fastq_list1_rand(fastq_list1, output_dir, rep_names_dict[rep])
        write_fastq_list2_rand(fastq_list2, output_dir, rep_names_dict[rep])
        counter = int(counter+(num_read/2 - int(read_count_dict[rep])))

#Write fastq_list1 random reads
def write_fastq_list1_rand(fastq_list1, output_dir, rep_name):
    read1 = os.path.join(output_dir, 'sim_' + str(rep_name) + '_read1.fastq')
    read1_write = open(read1, 'a')
    for index in range(len(fastq_list1)):
        for line in fastq_list1[index]:
            read1_write.write(line)
    read1_write.close()

#Write fastq_list2 random reads
def write_fastq_list2_rand(fastq_list2, output_dir, rep_name):
    read2 = os.path.join(output_dir, 'sim_' + str(rep_name) + '_read2.fastq')
    read2_write = open(read2, 'a')
    for index in range(len(fastq_list2)):
        for line in fastq_list2[index]:
            read2_write.write(line)
    read2_write.close()
