"""
Angel.py

This module creates the final report file.

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
import os.path
import statistics
import multiprocessing
from collections import OrderedDict
from subprocess import DEVNULL, STDOUT, check_call

#Create output info file
def output_info_total(intron_list, intron_dict, real_frag_dict, read_junc_dict, output_dir, strand_dict, fpkm_dict, intron_prop_dict, rep_names_dict):
    #Update some dicts
    for intron in intron_list:
        if intron[0] not in intron_prop_dict:
            intron_prop_dict[intron[0]] = 0.0
        if intron[0] not in fpkm_dict:
            fpkm_dict[intron[0]] = 0
    #Make intron_num_dict
    intron_num_dict = OrderedDict()
    for gene in intron_dict:
        for intron in intron_dict[gene]:
            try:
                intron_num_dict[gene].append(intron[6])
            except:
                intron_num_dict[gene]= []
                intron_num_dict[gene].append(intron[6])
    first_line = '#Gene_id\tTranscript_id\tIntron_no.\tCoordinates\tStrand\t'
    new_rep_names_dict = OrderedDict()
    for rep in rep_names_dict:
        new_rep_names_dict[rep] = rep_names_dict[rep][0].upper() + '_rep' + rep_names_dict[rep][-1]
    for rep in new_rep_names_dict:
        first_line = first_line + new_rep_names_dict[rep] + '_EE\t' + new_rep_names_dict[rep] + '_EI\t' + new_rep_names_dict[rep] + '_IE\t' + new_rep_names_dict[rep] + '_%IRratio\t'
    first_line = first_line + 'Simulated_ratio\tEstimated_gene_FPKM\n'
    output_path = os.path.join(output_dir, 'report.tsv')
    out = open(output_path, 'w')
    out.write(first_line)
    for intron in intron_list:
        outline = intron[0] + '\t' + intron[1] + '\t' + str(intron[6]) + '\t' + intron[2] + ':' + str(intron[3]) + '-' + str(intron[4]) + '\t' + strand_dict[intron[0]] + '\t'
        for rep in rep_names_dict:
            outline = outline + str(read_junc_dict[rep][str(intron[1]) + '-ee' + str(intron[6])]) + '\t' + str(read_junc_dict[rep][str(intron[1]) + '-ei' + str(intron[6])]) + '\t' + str(read_junc_dict[rep][str(intron[1]) + '-ie' + str(intron[6])]) + '\t' + str(round(statistics.mean([read_junc_dict[rep][str(intron[1]) + '-ei' + str(intron[6])], read_junc_dict[rep][str(intron[1]) + '-ie' + str(intron[6])]])/(statistics.mean([read_junc_dict[rep][str(intron[1]) + '-ei' + str(intron[6])], read_junc_dict[rep][str(intron[1]) + '-ie' + str(intron[6])]]) + read_junc_dict[rep][str(intron[1]) + '-ee' + str(intron[6])] + 0.000000001)*100, 1)) + '\t'
        if intron[0] in intron_num_dict:
            if intron[6] in intron_num_dict[intron[0]]:
                outline = outline + str(intron_prop_dict[intron[0]]) + '\t' + str(round(fpkm_dict[intron[0]], 3)) + '\n'
            else:
                outline = outline + '0.00\t' + str(round(fpkm_dict[intron[0]], 3)) + '\n'
        else:
            outline = outline + '0.00\t' + str(round(fpkm_dict[intron[0]], 3)) + '\n'
        out.write(outline)
    out.close()

#Pre-output info
def pre_output_info(intron_list, strand_dict, omit_gene_len_list):
    first_line = '#Gene_id\tTranscript_id\tIntron_no.\tCoordinates\tStrand\t'
    outline_dict = OrderedDict()
    for intron in intron_list:
        if intron[1] not in omit_gene_len_list:
            outline_dict[intron[0] + '-' + str(intron[6])] = intron[0] + '\t' + intron[1] + '\t' + str(intron[6]) + '\t' + intron[2] + ':' + str(intron[3]) + '-' + str(intron[4]) + '\t' + strand_dict[intron[0]] + '\t'
    return first_line, dict(outline_dict)

#Create output info per replicate
def output_info_rep(first_line, outline_dict, intron_list, read_junc_dict, rep_names_dict, omit_gene_len_list):
    new_rep_names_dict = OrderedDict()
    for rep in rep_names_dict:
        new_rep_names_dict[rep] = rep[0].upper() + '_rep' + rep_names_dict[rep][-1]
    for rep in new_rep_names_dict:
        first_line += new_rep_names_dict[rep] + '_EE\t' + new_rep_names_dict[rep] + '_EI\t' + new_rep_names_dict[rep] + '_IE\t' + new_rep_names_dict[rep] + '_%IRratio\t'
    for intron in intron_list:
        if intron[1] not in omit_gene_len_list:
            for rep in rep_names_dict:
                outline_dict[intron[0] + '-' + str(intron[6])] += str(read_junc_dict[rep][str(intron[1]) + '-ee' + str(intron[6])]) + '\t' + str(read_junc_dict[rep][str(intron[1]) + '-ei' + str(intron[6])]) + '\t' + str(read_junc_dict[rep][str(intron[1]) + '-ie' + str(intron[6])]) + '\t' + str(round(statistics.mean([read_junc_dict[rep][str(intron[1]) + '-ei' + str(intron[6])], read_junc_dict[rep][str(intron[1]) + '-ie' + str(intron[6])]])/(statistics.mean([read_junc_dict[rep][str(intron[1]) + '-ei' + str(intron[6])], read_junc_dict[rep][str(intron[1]) + '-ie' + str(intron[6])]]) + read_junc_dict[rep][str(intron[1]) + '-ee' + str(intron[6])] + 0.000000001)*100, 1)) + '\t'
    return first_line, dict(outline_dict)

#Create final output
def output_info_final(first_line, outline_dict, intron_list, intron_dict, output_dir, fpkm_dict, intron_prop_dict, intron_prop_dict_ctrl, omit_gene_len_list):
    #Update some dicts
    for intron in intron_list:
        if intron[0] not in intron_prop_dict:
            intron_prop_dict[intron[0]] = 0.0
            intron_prop_dict_ctrl[intron[0]] = 0.0
        if intron[0] not in fpkm_dict:
            fpkm_dict[intron[0]] = 0
    #Make intron_num_dict
    intron_num_dict = OrderedDict()
    for gene in intron_dict:
        for intron in intron_dict[gene]:
            try:
                intron_num_dict[gene].append(intron[6])
            except:
                intron_num_dict[gene]= []
                intron_num_dict[gene].append(intron[6])
    first_line += 'Simulated_ratio_A\tSimulated_ratio_B\tEstimated_gene_FPKM\n'
    output_path = os.path.join(output_dir, 'report.tsv')
    out = open(output_path, 'w')
    out.write(first_line)
    for intron in intron_list:
        if intron[1] not in omit_gene_len_list:
            if intron[0] in intron_num_dict:
                if intron[6] in intron_num_dict[intron[0]]:
                    outline_dict[intron[0] + '-' + str(intron[6])] += str(intron_prop_dict[intron[0]]) + '\t' + str(intron_prop_dict_ctrl[intron[0]]) + '\t' + str(round(fpkm_dict[intron[0]], 3)) + '\n'
                else:
                    outline_dict[intron[0] + '-' + str(intron[6])] += '0.00\t0.00\t' + str(round(fpkm_dict[intron[0]], 3)) + '\n'
            else:
                outline_dict[intron[0] + '-' + str(intron[6])] += '0.00\t0.00\t' + str(round(fpkm_dict[intron[0]], 3)) + '\n'
            out.write(outline_dict[intron[0] + '-' + str(intron[6])])
    out.close()

#Compress fastq multi
def compress_multi(rep_names_dict, output_dir, num_threads, total_reps):
    if int(num_threads) >= total_reps*2: #num_threads >= number of fastq files to be generated
        jobs = []
        for rep in rep_names_dict:
            path1 = os.path.join(output_dir, 'sim_' + str(rep_names_dict[rep]) + '_read1.fastq')
            path2 = os.path.join(output_dir, 'sim_' + str(rep_names_dict[rep]) + '_read2.fastq')
            p = multiprocessing.Process(target=make_gz, args=(path1,))
            jobs.append(p)
            p.start()
            p = multiprocessing.Process(target=make_gz, args=(path2,))
            jobs.append(p)
            p.start()
        for job in jobs:
            job.join()
        jobs = []
        for rep in rep_names_dict:
            path1 = os.path.join(output_dir, 'sim_' + str(rep_names_dict[rep]) + '_read1.fastq')
            path2 = os.path.join(output_dir, 'sim_' + str(rep_names_dict[rep]) + '_read2.fastq')
            p = multiprocessing.Process(target=delete_fastq, args=(path1,))
            jobs.append(p)
            p.start()
            p = multiprocessing.Process(target=delete_fastq, args=(path2,))
            jobs.append(p)
            p.start()
        for job in jobs:
            job.join()
    elif int(num_threads) < total_reps*2:
        for rep in rep_names_dict:
            path1 = os.path.join(output_dir, 'sim_' + str(rep_names_dict[rep]) + '_read1.fastq')
            path2 = os.path.join(output_dir, 'sim_' + str(rep_names_dict[rep]) + '_read2.fastq')
            make_gz(path1)
            make_gz(path2)
            delete_fastq(path1)
            delete_fastq(path2)

#Compress gzip file
def make_gz(path):
    f = open(path, 'rb')
    fastq = f.read()
    f.close()
    out = gzip.open(path + '.gz', 'wb')
    out.write(fastq)
    out.close()

#Delete fastq files
def delete_fastq(path):
    check_call(['rm', path], stdout=DEVNULL, stderr=STDOUT)
