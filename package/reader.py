"""
Reader.py

This module extracts FASTQ reads from DWGSIM.

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
import numpy as np
import multiprocessing
from collections import OrderedDict

#For sorting by read id
def take_read_id(elem):
    return elem[0].rsplit('_', 1)[1].split('/')[0]

#Extract fastq and sort strandness (Multi_thread, Pickle)
def extract_fastq_pickle(fastq_name, output_dir, trans_gene_dict, strand_dict):
    fastq1 = os.path.join(output_dir, fastq_name + '.read1.fastq.gz')
    fastq2 = os.path.join(output_dir, fastq_name + '.read2.fastq.gz')
    fastq_dict1 = OrderedDict()
    fastq_dict2 = OrderedDict()
    read_name_list = []
    with gzip.open(fastq1, 'rt') as f1:
        for line in f1:
            if line.rsplit('_', 9)[0] not in read_name_list:
                read_name_list.append(line.rsplit('_', 9)[0])
            #Sort reads according to strandedness
            if strand_dict[trans_gene_dict[line.rsplit('_', 9)[0][1:-2]]] == '+':
                #If read1 is forward
                if line.rsplit('_', 7)[1] == '0':
                    #Place in dict1
                    try:
                        #dum = fastq_dict1[line.rsplit('_', 9)[0]]
                        fastq_dict1[line.rsplit('_', 9)[0]].append([line, next(f1), next(f1), next(f1)])
                    except:
                        fastq_dict1[line.rsplit('_', 9)[0]] = []
                        fastq_dict1[line.rsplit('_', 9)[0]].append([line, next(f1), next(f1), next(f1)])
                #Elif read1 is reverse
                elif line.rsplit('_', 7)[1] == '1':
                    #Place in dict2 and change read number
                    try:
                        #dum = fastq_dict2[line.rsplit('_', 9)[0]]
                        fastq_dict2[line.rsplit('_', 9)[0]].append([line.rsplit('/', 1)[0] + '/2\n', next(f1), next(f1), next(f1)])
                    except:
                        fastq_dict2[line.rsplit('_', 9)[0]] = []
                        fastq_dict2[line.rsplit('_', 9)[0]].append([line.rsplit('/', 1)[0] + '/2\n', next(f1), next(f1), next(f1)])
            elif strand_dict[trans_gene_dict[line.rsplit('_', 9)[0][1:-2]]] == '-':
                #If read1 is forward
                if line.rsplit('_', 7)[1] == '0':
                    #Place in dict2 and change read number
                    try:
                        #dum = fastq_dict2[line.rsplit('_', 9)[0]]
                        fastq_dict2[line.rsplit('_', 9)[0]].append([line.rsplit('/', 1)[0] + '/2\n', next(f1), next(f1), next(f1)])
                    except:
                        fastq_dict2[line.rsplit('_', 9)[0]] = []
                        fastq_dict2[line.rsplit('_', 9)[0]].append([line.rsplit('/', 1)[0] + '/2\n', next(f1), next(f1), next(f1)])
                #Elif read1 is reverse
                elif line.rsplit('_', 7)[1] == '1':
                    #Place in dict1
                    try:
                        #dum = fastq_dict1[line.rsplit('_', 9)[0]]
                        fastq_dict1[line.rsplit('_', 9)[0]].append([line, next(f1), next(f1), next(f1)])
                    except:
                        fastq_dict1[line.rsplit('_', 9)[0]] = []
                        fastq_dict1[line.rsplit('_', 9)[0]].append([line, next(f1), next(f1), next(f1)])
    with gzip.open(fastq2, 'rt') as f2:
        for line in f2:
            #Sort reads according to strandedness
            if strand_dict[trans_gene_dict[line.rsplit('_', 9)[0][1:-2]]] == '+':
                #If read1 is forward
                if line.rsplit('_', 7)[1] == '0':
                    #Place in dict2
                    try:
                        #dum = fastq_dict2[line.rsplit('_', 9)[0]]
                        fastq_dict2[line.rsplit('_', 9)[0]].append([line, next(f2), next(f2), next(f2)])
                    except:
                        fastq_dict2[line.rsplit('_', 9)[0]] = []
                        fastq_dict2[line.rsplit('_', 9)[0]].append([line, next(f2), next(f2), next(f2)])
                #Elif read1 is reverse
                elif line.rsplit('_', 7)[1] == '1':
                    #Place in dict1 and change read number
                    try:
                        #dum = fastq_dict1[line.rsplit('_', 9)[0]]
                        fastq_dict1[line.rsplit('_', 9)[0]].append([line.rsplit('/', 1)[0] + '/1\n', next(f2), next(f2), next(f2)])
                    except:
                        fastq_dict1[line.rsplit('_', 9)[0]] = []
                        fastq_dict1[line.rsplit('_', 9)[0]].append([line.rsplit('/', 1)[0] + '/1\n', next(f2), next(f2), next(f2)])
            elif strand_dict[trans_gene_dict[line.rsplit('_', 9)[0][1:-2]]] == '-':
                #If read1 is forward
                if line.rsplit('_', 7)[1] == '0':
                    #Place in dict1 and change read number
                    try:
                        #dum = fastq_dict1[line.rsplit('_', 9)[0]]
                        fastq_dict1[line.rsplit('_', 9)[0]].append([line.rsplit('/', 1)[0] + '/1\n', next(f2), next(f2), next(f2)])
                    except:
                        fastq_dict1[line.rsplit('_', 9)[0]] = []
                        fastq_dict1[line.rsplit('_', 9)[0]].append([line.rsplit('/', 1)[0] + '/1\n', next(f2), next(f2), next(f2)])
                #Elif read1 is reverse
                elif line.rsplit('_', 7)[1] == '1':
                    #Place in dict2
                    try:
                        #dum = fastq_dict2[line.rsplit('_', 9)[0]]
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

#Multiprocessing extract_fastq (Pickle) stranded
def multi_thread_extract_fastq_pickle(current_gene_list, output_dir, trans_gene_dict, strand_dict, num_threads):    
    #manager = multiprocessing.Manager()
    #dict1 = manager.dict()
    #dict2 = manager.dict()
    jobs = []
    if len(current_gene_list) >= num_threads:
        for index in range(int(num_threads)):
            p = multiprocessing.Process(target=extract_fastq_pickle, args=('tmp' + str(index) + '.fa.out.bwa', output_dir, trans_gene_dict, strand_dict))
            jobs.append(p)
            p.start()
    else:
        for index in range(len(current_gene_list)):
            p = multiprocessing.Process(target=extract_fastq_pickle, args=('tmp' + str(index) + '.fa.out.bwa', output_dir, trans_gene_dict, strand_dict))
            jobs.append(p)
            p.start()
    for job in jobs:
        job.join()
    fastq_dict1 = OrderedDict()
    fastq_dict2 = OrderedDict()
    #Unpickle
    if len(current_gene_list) >= num_threads:
        for index in range(int(num_threads)):
            pickle1 = open(os.path.join(output_dir, 'tmp' + str(index) + '.fa.out.bwa') + '1.pkl', 'rb')
            pickle2 = open(os.path.join(output_dir, 'tmp' + str(index) + '.fa.out.bwa') + '2.pkl', 'rb')
            dict1 = pickle.load(pickle1)
            dict2 = pickle.load(pickle2)
            pickle1.close()
            pickle2.close()
            fastq_dict1.update(dict1)
            fastq_dict2.update(dict2)
    else:
        for index in range(len(current_gene_list)):
            pickle1 = open(os.path.join(output_dir, 'tmp' + str(index) + '.fa.out.bwa') + '1.pkl', 'rb')
            pickle2 = open(os.path.join(output_dir, 'tmp' + str(index) + '.fa.out.bwa') + '2.pkl', 'rb')
            dict1 = pickle.load(pickle1)
            dict2 = pickle.load(pickle2)
            pickle1.close()
            pickle2.close()
            fastq_dict1.update(dict1)
            fastq_dict2.update(dict2)
    return dict(fastq_dict1), dict(fastq_dict2)

#Extract fastq and do not sort strandness (Multi_thread, Pickle)
def extract_fastq_pickle_not_stranded(fastq_name, output_dir, trans_gene_dict):
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

#Multiprocessing extract_fastq (Pickle) not_stranded
def multi_thread_extract_fastq_pickle_not_stranded(output_dir, trans_gene_dict, strand_dict, num_threads):    
    #manager = multiprocessing.Manager()
    #dict1 = manager.dict()
    #dict2 = manager.dict()
    jobs = []
    for index in range(int(num_threads)):
        p = multiprocessing.Process(target=extract_fastq_pickle_not_stranded, args=('tmp' + str(index) + '.fa.out.bwa', output_dir, trans_gene_dict))
        jobs.append(p)
        p.start()
    for job in jobs:
        job.join()
    fastq_dict1 = OrderedDict()
    fastq_dict2 = OrderedDict()
    #Unpickle
    for index in range(int(num_threads)):
        pickle1 = open(os.path.join(output_dir, 'tmp' + str(index) + '.fa.out.bwa') + '1.pkl', 'rb')
        pickle2 = open(os.path.join(output_dir, 'tmp' + str(index) + '.fa.out.bwa') + '2.pkl', 'rb')
        dict1 = pickle.load(pickle1)
        dict2 = pickle.load(pickle2)
        pickle1.close()
        pickle2.close()
        fastq_dict1.update(dict1)
        fastq_dict2.update(dict2)
    return dict(fastq_dict1), dict(fastq_dict2)

#Make fragment_dict
def make_fragment_dict(fastq_dict, intron_dict, intron_prop_dict, fpkm_dict, trans_gene_dict, num_read, cdna_size_dict, max_rep_var, dataset_type, ran_read, intron_size_dict):
    transcript_frag_dict = OrderedDict()
    real_frag_dict = OrderedDict()
    #Create dictionary of percent variance for each transcript
    transcript_perc_var_dict = OrderedDict()
    for transcript_id in fastq_dict:
        transcript = transcript_id[1:].split('-')[0]
        transcript_perc_var_dict[transcript] = np.random.uniform(-max_rep_var/100,max_rep_var/100) + 1
    if dataset_type == 'intron':
        for transcript_id in fastq_dict:
            transcript = transcript_id[1:].split('-')[0]
            #Check if transcript is exonic or intronic
            if transcript_id[-1] == 'e':
                #Check if transcript has intron partner
                if trans_gene_dict[transcript] in intron_dict:
                    transcript_frag_dict[transcript_id] = round((((fpkm_dict[trans_gene_dict[transcript]]*(cdna_size_dict[transcript]/1000))*(num_read*(1-ran_read)/2/1000000))*(1-intron_prop_dict[trans_gene_dict[transcript]]))*transcript_perc_var_dict[transcript])
                    real_frag_dict[transcript_id] = round((((fpkm_dict[trans_gene_dict[transcript]]*(cdna_size_dict[transcript]/1000))*(num_read*(1-ran_read)/2/1000000))*(1-intron_prop_dict[trans_gene_dict[transcript]]))*transcript_perc_var_dict[transcript])
                else:
                    transcript_frag_dict[transcript_id] = round(((fpkm_dict[trans_gene_dict[transcript]]*(cdna_size_dict[transcript]/1000))*(num_read*(1-ran_read)/2/1000000))*transcript_perc_var_dict[transcript])
                    real_frag_dict[transcript_id] = round(((fpkm_dict[trans_gene_dict[transcript]]*(cdna_size_dict[transcript]/1000))*(num_read*(1-ran_read)/2/1000000))*transcript_perc_var_dict[transcript])
            elif transcript_id[-1] == 'i':
                #Decipher number of fragments
                transcript_frag_dict[transcript_id] = round((((fpkm_dict[trans_gene_dict[transcript]]*(cdna_size_dict[transcript]/1000))*(num_read*(1-ran_read)/2/1000000))*intron_prop_dict[trans_gene_dict[transcript]])*transcript_perc_var_dict[transcript]*(intron_size_dict[transcript]/cdna_size_dict[transcript]))
                real_frag_dict[transcript_id] = round((((fpkm_dict[trans_gene_dict[transcript]]*(cdna_size_dict[transcript]/1000))*(num_read*(1-ran_read)/2/1000000))*intron_prop_dict[trans_gene_dict[transcript]])*transcript_perc_var_dict[transcript])
    elif dataset_type == 'control':
        for transcript_id in fastq_dict:
            transcript = transcript_id[1:].split('-')[0]
            #Check if transcript is exonic or intronic
            if transcript_id[-1] == 'e':
                transcript_frag_dict[transcript_id] = round(((fpkm_dict[trans_gene_dict[transcript]]*(cdna_size_dict[transcript]/1000))*(num_read*(1-ran_read)/1000000)/2)*transcript_perc_var_dict[transcript])
                real_frag_dict[transcript_id] = round(((fpkm_dict[trans_gene_dict[transcript]]*(cdna_size_dict[transcript]/1000))*(num_read*(1-ran_read)/1000000)/2)*transcript_perc_var_dict[transcript])
    return dict(transcript_frag_dict), dict(real_frag_dict)
