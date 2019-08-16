"""

Compose.py

This module writes FASTQ reads.

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

import pickle
import os.path
import multiprocessing
from collections import OrderedDict

#Make output files
def make_out_files_many(output_dir, transcript_frag_dict, intron_reps, control_reps):
    rep_names = []
    rep_names_dict = OrderedDict()
    for num in range(intron_reps):
        rep_names.append('intron_rep' + str(num + 1))
    for num in range(control_reps):
        rep_names.append('control_rep' + str(num + 1))    
    for rep_name in rep_names:
        read1 = os.path.join(output_dir, 'sim_' + rep_name + '_read1.fastq')
        read2 = os.path.join(output_dir, 'sim_' + rep_name + '_read2.fastq')
        read1_write = open(read1, 'w')
        read2_write = open(read2, 'w')
        read1_write.close()
        read2_write.close()
    for rep in transcript_frag_dict:
        rep_names_dict[rep] = rep_names[int(rep.split('-')[1])]
    return dict(rep_names_dict)

#Make output files
def make_out_files(output_dir, transcript_frag_dict, replicate):
    if replicate[0] == 'A':
        rep_name = 'sample_A_rep' + str(int(replicate[-1]) + 1)
    elif replicate[0] == 'B':
        rep_name = 'sample_B_rep' + str(int(replicate[-1]) + 1)
    read1 = os.path.join(output_dir, 'sim_' + rep_name + '_read1.fastq')
    read2 = os.path.join(output_dir, 'sim_' + rep_name + '_read2.fastq')
    read1_write = open(read1, 'w')
    read2_write = open(read2, 'w')
    read1_write.close()
    read2_write.close()
    rep_names_dict = OrderedDict()
    for rep in transcript_frag_dict:
        rep_names_dict[rep] = rep_name
    return dict(rep_names_dict)

#Make read index chunks for each transcript
def read_index_chunks(transcript_frag_dict, fastq_dict1):
    chunks_dict = OrderedDict()
    latest_index = OrderedDict()
    index = 0
    for rep in transcript_frag_dict:
        chunks_dict[rep] = OrderedDict()
    for transcript_id in fastq_dict1:
        latest_index[transcript_id] = 0 
    for rep in transcript_frag_dict:
        for transcript_id in transcript_frag_dict[rep]:
            #try: #At low transcript simulation amounts, sometimes DWGSIM does not simulate some transcripts.
            chunks_dict[rep][transcript_id] = [latest_index[transcript_id], min(transcript_frag_dict[rep][transcript_id] + latest_index[transcript_id], len(fastq_dict1[transcript_id]))]
            latest_index[transcript_id] = min(transcript_frag_dict[rep][transcript_id] + latest_index[transcript_id], len(fastq_dict1[transcript_id]))
            #except:
                #chunks_dict[rep][transcript_id] = [0, 0]
            #except:
                #chunks_dict[rep][transcript_id] = []
                #chunks_dict[rep][transcript_id].append([latest_index[transcript_id], min(transcript_frag_dict[rep][transcript_id] + latest_index[transcript_id], len(fastq_dict1[transcript_id]))])
    return dict(chunks_dict)

#Update Fastq lists
def update_fastq(fastq_list1, fastq_list2, transcript_frag_dict, fastq_dict1, fastq_dict2, read_index_chunks):
    remove_list = []
    for transcript_id in transcript_frag_dict:
        start = read_index_chunks[transcript_id][0]
        end = read_index_chunks[transcript_id][1]
        #Split fastq_dict[transcript_id] index into chunks dict
        fastq_list1.extend(fastq_dict1[transcript_id][start:end])
        fastq_list2.extend(fastq_dict2[transcript_id][start:end])
        length = len(fastq_dict1[transcript_id][start:end])
        #Check if there were enough reads
        if transcript_frag_dict[transcript_id] == length:
            remove_list.append(transcript_id)
        elif transcript_frag_dict[transcript_id] > length:
            transcript_frag_dict[transcript_id] = transcript_frag_dict[transcript_id] - length
        else:
            print("ERROR: Update_fastq error")
            break
    for transcript_id in remove_list:
        del transcript_frag_dict[transcript_id]

#Write fastq_list1
def write_fastq_list1(fastq_list1, output_dir, rep_name, junc_dict, total_intron_dict, read_length):
    read1 = os.path.join(output_dir, 'sim_' + str(rep_name) + '_read1.fastq')
    read1_write = open(read1, 'a')
    temp_range_dict = OrderedDict()
    for index in range(len(fastq_list1)):
        n = 1
        for line in fastq_list1[index]:
            read1_write.write(line)
            if n == 1:
                try:
                    temp_range_dict[line.rsplit('_', 9)[0]].append([int(line.rsplit('_', 9)[1]), int(line.rsplit('_', 9)[1]) + int(read_length)]) #count_junc(line, str(rep_name), junc_dict)
                except:
                    temp_range_dict[line.rsplit('_', 9)[0]] = []
                    temp_range_dict[line.rsplit('_', 9)[0]].append([int(line.rsplit('_', 9)[1]), int(line.rsplit('_', 9)[1]) + int(read_length)])
            n += 1
        n = 1
    read1_write.close()
    temp_junc_dict = count_junc(temp_range_dict, rep_name, junc_dict, total_intron_dict)
    pickle1 = open(os.path.join(output_dir, rep_name) + '1.pkl', 'wb')
    pickle.dump(temp_junc_dict, pickle1, -1)
    pickle1.close()

#Write fastq_list2
def write_fastq_list2(fastq_list2, output_dir, rep_name, junc_dict, total_intron_dict, read_length):
    read2 = os.path.join(output_dir, 'sim_' + str(rep_name) + '_read2.fastq')
    read2_write = open(read2, 'a')
    temp_range_dict = OrderedDict()
    for index in range(len(fastq_list2)):
        n = 1
        for line in fastq_list2[index]:
            read2_write.write(line)
            if n == 1:
                try:
                    temp_range_dict[line.rsplit('_', 9)[0]].append([int(line.rsplit('_', 9)[2]), int(line.rsplit('_', 9)[2]) + int(read_length)]) #count_junc(line, str(rep_name), junc_dict)
                except:
                    temp_range_dict[line.rsplit('_', 9)[0]] = []
                    temp_range_dict[line.rsplit('_', 9)[0]].append([int(line.rsplit('_', 9)[2]), int(line.rsplit('_', 9)[2]) + int(read_length)])
            n += 1
        n = 1
    read2_write.close()
    temp_junc_dict = count_junc(temp_range_dict, rep_name, junc_dict, total_intron_dict)
    pickle2 = open(os.path.join(output_dir, rep_name) + '2.pkl', 'wb')
    pickle.dump(temp_junc_dict, pickle2, -1)
    pickle2.close()

#Count reads at splicing junctions
def count_junc(temp_range_dict, rep_name, junc_dict, total_intron_dict):
    #Set overlap threshold
    n = 8
    temp_junc_dict = OrderedDict()
    for transcript in total_intron_dict:
        for num in total_intron_dict[transcript]:
            temp_junc_dict[transcript + '-ee' + str(num)] = 0
            temp_junc_dict[transcript + '-ei' + str(num)] = 0
            temp_junc_dict[transcript + '-ie' + str(num)] = 0
    for transcript_id in temp_range_dict:
        if transcript_id[1:-2] in total_intron_dict: #If transcript has an intron
            if transcript_id[-1] == 'e':
                for num in total_intron_dict[transcript_id[1:-2]]:
                    for rrange in temp_range_dict[transcript_id]:
                        if rrange[0]+n <= junc_dict[transcript_id[1:-2] + '-ee' + str(num)] <= rrange[1]-n:
                            temp_junc_dict[transcript_id[1:-2] + '-ee' + str(num)] += 1
            elif transcript_id[-1] == 'i':
                for num in total_intron_dict[transcript_id[1:-2]]:
                    for rrange in temp_range_dict[transcript_id]:
                        if rrange[0]+n <= junc_dict[transcript_id[1:-2] + '-ei' + str(num)] <= rrange[1]-n:
                            temp_junc_dict[transcript_id[1:-2] + '-ei' + str(num)] += 1
                        if rrange[0]+n <= junc_dict[transcript_id[1:-2] + '-ie' + str(num)] <= rrange[1]-n:
                            temp_junc_dict[transcript_id[1:-2] + '-ie' + str(num)] += 1
    return dict(temp_junc_dict)

#Prepare objects for count_junc
def junc_count_prep(intron_list, junc_dict, rep_names_dict):
    read_junc_dict = OrderedDict()
    for rep in rep_names_dict:
        read_junc_dict[rep] = OrderedDict()
        for intron in intron_list:
            read_junc_dict[rep][str(intron[1]) + '-ee' + str(intron[6])] = 0
            read_junc_dict[rep][str(intron[1]) + '-ei' + str(intron[6])] = 0
            read_junc_dict[rep][str(intron[1]) + '-ie' + str(intron[6])] = 0
    total_intron_dict = OrderedDict()
    for intron in intron_list:
        junc_dict[intron[1] + '-ee' + str(intron[6])] = intron[5]
        if str(intron[1] + '-ei' + str(intron[6])) not in junc_dict:
            junc_dict[intron[1] + '-ei' + str(intron[6])] = -1
            junc_dict[intron[1] + '-ie' + str(intron[6])] = -1
        try:
            total_intron_dict[intron[1]].append(intron[6])
        except:
            total_intron_dict[intron[1]] = []
            total_intron_dict[intron[1]].append(intron[6])
    return dict(read_junc_dict), dict(total_intron_dict)

#Multi-process fastq_list and write fastq
def multi_fastq_list_write(num_threads, transcript_frag_dict, fastq_dict1, fastq_dict2, total_reps, output_dir, rep_names_dict, chunks_dict, junc_dict, read_junc_dict, total_intron_dict, read_length):
    #Create fastq lists
    fastq_list1 = OrderedDict()
    fastq_list2 = OrderedDict()
    if int(num_threads) >= total_reps*2: #num_threads >= number of fastq files to be generated
        jobs = []
        for rep in transcript_frag_dict:
            fastq_list1[rep] = []
            fastq_list2[rep] = []
            update_fastq(fastq_list1[rep], fastq_list2[rep], transcript_frag_dict[rep], fastq_dict1, fastq_dict2, chunks_dict[rep])
            p = multiprocessing.Process(target=write_fastq_list1, args=(fastq_list1[rep], output_dir, rep_names_dict[rep], junc_dict, total_intron_dict, read_length))
            jobs.append(p)
            p.start()
            p = multiprocessing.Process(target=write_fastq_list2, args=(fastq_list2[rep], output_dir, rep_names_dict[rep], junc_dict, total_intron_dict, read_length))
            jobs.append(p)
            p.start()
        for job in jobs:
            job.join()
    elif int(num_threads) < total_reps*2:
        #No multiprocessing
        for rep in transcript_frag_dict:
            fastq_list1[rep] = []
            fastq_list2[rep] = []
            update_fastq(fastq_list1[rep], fastq_list2[rep], transcript_frag_dict[rep], fastq_dict1, fastq_dict2, chunks_dict[rep])
            write_fastq_list1(fastq_list1[rep], output_dir, rep_names_dict[rep], junc_dict, total_intron_dict)
            write_fastq_list2(fastq_list2[rep], output_dir, rep_names_dict[rep], junc_dict, total_intron_dict)
    for rep in transcript_frag_dict:
        pickle1 = open(os.path.join(output_dir, rep_names_dict[rep]) + '1.pkl', 'rb')
        pickle2 = open(os.path.join(output_dir, rep_names_dict[rep]) + '2.pkl', 'rb')
        dict1 = pickle.load(pickle1)
        dict2 = pickle.load(pickle2)
        pickle1.close()
        pickle2.close()
        for i in dict1:
            read_junc_dict[rep][i] += dict1[i]
        for i in dict2:
            read_junc_dict[rep][i] += dict2[i]

#Update fasta_dict
def update_fasta_dict(fasta_dict, transcript_frag_dict, trans_gene_dict):
    new_len = 0
    new_gene_list = []
    remove_list = []
    #Get total required transcipts
    total_transcripts = []
    for rep in transcript_frag_dict:
        for transcript_id in transcript_frag_dict[rep]:
            if transcript_id not in total_transcripts:
                total_transcripts.append(transcript_id)
    for transcript in fasta_dict:
        if "@"+transcript not in total_transcripts:
            remove_list.append(transcript)
        else:
            new_len = new_len + len(fasta_dict[transcript])
            new_gene_list.append(trans_gene_dict[transcript[:-2]])
    for transcript in remove_list:
        del fasta_dict[transcript]
    return new_len, new_gene_list
