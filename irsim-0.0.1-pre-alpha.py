#!/usr/bin/python3

"""
IRSim: Simulation of Intron Retention in Coding RNA

Version 0.0.1-pre-alpha

This is the main script of the program IRSim.

Copyright (C) 2019 Tham Cheng Yong

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import sys
import math
import copy
import time
import gzip
import pickle
import os.path
import datetime
import statistics 
import numpy as np
import configparser
from sys import argv
import multiprocessing
from pyfaidx import Fasta
from collections import OrderedDict
from subprocess import DEVNULL, STDOUT, check_call

#Check sections and keys of config.ini file
def config_checker(config):
    for key in ['Output directory', 'DWGSIM directory']:
        if config.has_option('Paths', key) == False:
            print('Error in config.ini: "%s" key do not exists' % (key))
            sys.exit()
    for key in ['Number of replicates for dataset with intron retention', 'Number of replicates for dataset without intron retention', 'Maximum percent variance between replicates', 'Number of reads per replicate', 'Paired-end sequencing', 'Strand-specific', 'Read length', 'Average outer distance between read pairs (insert length)', 'Standard deviation of outer distance between read pairs', 'Per base error rate of the first read', 'Per base error rate of the second read', 'Rate of mutation', 'Fraction of mutations that are indels', 'Probability an indel is extended', 'Minimum length indel', 'Probability of a random read']:
        if config.has_option('Sequencing-details', key) == False:
            print('Error in config.ini: "%s" key do not exists' % (key))
            sys.exit()
    for key in ['Percentage of total introns which are retained', 'Probability of retaining an intron with length more than 1000 bases', 'Probability of retaining an intron with length between 100-1000 bases', 'Probability of retaining an intron with length less than 100 bases', 'Proportions of transcripts with intron retention']:
        if config.has_option('Intron-retention', key) == False:
            print('Error in config.ini: "%s" key do not exists' % (key))
            sys.exit()
    for key in ['Lower limit of FPKM in gene expression model', 'Upper limit of FPKM in gene expression model']:
        if config.has_option('Gene-expression', key) == False:
            print('Error in config.ini: "%s" key do not exists' % (key))
            sys.exit()
    for key in ['Random seed']:
        if config.has_option('Seed', key) == False:
            print('Error in config.ini: "%s" key do not exists' % (key))
            sys.exit()
    for key in ['Number of threads']:
        if config.has_option('Threads', key) == False:
            print('Error in config.ini: "%s" key do not exists' % (key))
            sys.exit()

#Extract configured seed value
def assign_seed(seed):
    try:
        np.random.seed(int(seed))
    except:
        np.random.seed(None)

#Check output dir exist
def out_dir_exist(output_dir):
    if not os.path.exists(output_dir):
        sys.exit("ERROR: Output directory '%s' is not found." % (output_dir))

#Extract FPKM_model gene expression
def ge_extract(model_file, minimum_FPKM, maximum_FPKM):
    fpkm_dict = OrderedDict()
    with open(model_file) as model:
        lowest_FPKM = maximum_FPKM
        highest_FPKM = minimum_FPKM
        for line in model:
            #Check if FPKM is within boundaries
            if minimum_FPKM <= float(line.split('\t')[1]) <= maximum_FPKM:
                #Make gene dictionary
                fpkm_dict[line.split('\t')[0]] = float(line.split('\t')[1])
            #Get lowest_FPKM
            if float(line.split('\t')[1]) != 0 and float(line.split('\t')[1]) >= minimum_FPKM and float(line.split('\t')[1]) < lowest_FPKM:
                lowest_FPKM = float(line.split('\t')[1])
            #Get highest_FPKM
            if float(line.split('\t')[1]) <= maximum_FPKM and float(line.split('\t')[1]) > highest_FPKM:
                highest_FPKM = float(line.split('\t')[1])
    if minimum_FPKM < 0:
        sys.exit("ERROR: Assigned lower limit FPKM value %s is less than 0." % (minimum_FPKM))
    return dict(fpkm_dict), lowest_FPKM, highest_FPKM

#Obtain intron information (The first transcript of each gene ordered in the GTF file is selected as the only transcript)
def get_intron_info(gtf_file, fpkm_dict):
    gtf = open(gtf_file, 'r').read().splitlines()
    num_lines = len(gtf)
    #Search for first entry and append it to last line in gtf as dummy
    for line in gtf:
        #Ignore comment lines
        if not line.strip().startswith("#"):
            gtf.append(line)
            break
    exon_tmp_list = []
    gene_list = []
    trans_gene_dict = OrderedDict()
    transcript_dict = OrderedDict()
    strand_dict = OrderedDict()
    biotype_dict = OrderedDict()
    biotype_dict_transcript = OrderedDict()
    intron_list = []
    intron_size = 0
    transcript_coord_dict = OrderedDict()
    for index in range(num_lines):
        #Ignore comment lines
        if not gtf[index].strip().startswith("#"):
            #Check next line for same gene_id and make exon list(Ensure GTF uses double quotes and correct columns) Col: gene_id,transcript_id,chr,start,end,strand
            if gtf[index].split('\t')[8].split(';')[0].split('"')[1] == gtf[index+1].split('\t')[8].split(';')[0].split('"')[1]:
                if gtf[index].split('\t')[2] == "exon":
                    exon_tmp_list.append([str(gtf[index].split('\t')[8].split(';')[0].split('"')[1]), str(gtf[index].split('\t')[8].split(';')[1].split('"')[1]), str(gtf[index].split('\t')[0]), int(gtf[index].split('\t')[3]), int(gtf[index].split('\t')[4]), str(gtf[index].split('\t')[6])])
                    biotype_dict[gtf[index].split('\t')[8].split(';')[0].split('"')[1]] = gtf[index].split('\t')[8].split(';')[5].split('"')[1]
            else:
                if gtf[index].split('\t')[2] == "exon":
                    exon_tmp_list.append([str(gtf[index].split('\t')[8].split(';')[0].split('"')[1]), str(gtf[index].split('\t')[8].split(';')[1].split('"')[1]), str(gtf[index].split('\t')[0]), int(gtf[index].split('\t')[3]), int(gtf[index].split('\t')[4]), str(gtf[index].split('\t')[6])])
                    biotype_dict[gtf[index].split('\t')[8].split(';')[0].split('"')[1]] = gtf[index].split('\t')[8].split(';')[5].split('"')[1]
                #Check if gene is in selected list
                if exon_tmp_list[0][0] in fpkm_dict:
                    #Check if gene is protein coding
                    if biotype_dict[exon_tmp_list[0][0]] == 'protein_coding':
                        #Append new gene list
                        gene_list.append(exon_tmp_list[0][0])
                        #Choose first transcript_id and discard the rest
                        first_transcript = exon_tmp_list[0][1]
                        remove_list = []
                        for i in exon_tmp_list:
                            if i[1] != first_transcript:
                                remove_list.append(i)
                        for i in remove_list:
                            exon_tmp_list.remove(i)
                        #Make gene transcript dict
                        trans_gene_dict[first_transcript] = exon_tmp_list[0][0]
                        #Make gene transcript dict
                        transcript_dict[exon_tmp_list[0][0]] = first_transcript
                        #Make gene strand dict
                        strand_dict[exon_tmp_list[0][0]] = exon_tmp_list[0][5]
                        #Obtain intron information Col: gene_id,transcript_id,chr,start,end,transcript_ins_point,intron_num,intron_size
                        #Check strandness
                        if exon_tmp_list[0][5] == '+':
                            tsc = exon_tmp_list[0][3] #transcript_start_coord
                            exon_num = len(exon_tmp_list)
                            for index in range(exon_num-1):
                                intron_list.append([exon_tmp_list[index][0], exon_tmp_list[index][1], exon_tmp_list[index][2], int(exon_tmp_list[index][4]+1), int(exon_tmp_list[index+1][3]-1), int(exon_tmp_list[index][4]+1-tsc-intron_size), index+1, int(exon_tmp_list[index+1][3]-1-(exon_tmp_list[index][4]+1)+1)])
                                intron_size = intron_size + int(exon_tmp_list[index+1][3]-1-(exon_tmp_list[index][4]+1)+1)
                        elif exon_tmp_list[0][5] == '-':
                            tsc = exon_tmp_list[-1][3] #transcript_start_coord
                            exon_num = len(exon_tmp_list)
                            for index in reversed(range(1,exon_num)):
                                intron_list.append([exon_tmp_list[index][0], exon_tmp_list[index][1], exon_tmp_list[index][2], int(exon_tmp_list[index][4]+1), int(exon_tmp_list[index-1][3]-1), int(exon_tmp_list[index][4]+1-tsc-intron_size), index, int(exon_tmp_list[index-1][3]-1-(exon_tmp_list[index][4]+1)+1)])
                                intron_size = intron_size + int(exon_tmp_list[index-1][3]-1-(exon_tmp_list[index][4]+1)+1)
                        else:
                            sys.exit("ERROR: Gene %s strandness neither + nor -" % (exon_tmp_list[0][0]))
                intron_size = 0
                exon_tmp_list = []
    #Get transcript coordinates
    for line in gtf:
        if not line.strip().startswith("#"):
            if line.split('\t')[2] == 'transcript':
                transcript_coord_dict[str(line.split('\t')[8].split(';')[1].split('"')[1])] = line.split('\t')[0] + ':' + line.split('\t')[3] + '-' + line.split('\t')[4] + ';' + line.split('\t')[6]
    return gene_list, intron_list, dict(transcript_dict), dict(strand_dict), dict(biotype_dict), dict(trans_gene_dict), dict(transcript_coord_dict)

#Take seventh element
def takeSeventh(elem):
    return int(elem[6])

#Reverse complement
tab = str.maketrans("ACTG","TGAC")
def rc(seq):
    return seq.upper().translate(tab)[::-1]

#Select introns to retain based on intron sizes probabilities (Ptop: >1000 bases, mid: 100-1000 bases, bot: <100 bases intron sizes)
def select_introns(intron_list, intron_perc, top, mid, bot, prop_list, strand_dict):
    selected_introns = []
    intron_prop_dict = OrderedDict()
    num_introns = len(intron_list)
    num_cap = round(num_introns*intron_perc)
    num = 0
    while num < num_cap:
        np.random.shuffle(intron_list)
        for i in intron_list:
            if i not in selected_introns:
                if i[-1] > 1000:
                    if np.random.uniform(0.0, 1.0) >= top: #***bug*** it should be <= (Fixed in 0.0.2)
                        selected_introns.append(i)
                        intron_prop_dict[i[0]] = float(np.random.choice(prop_list).strip())/100
                        num += 1
                elif 100 <= i[-1] <= 1000:
                    if np.random.uniform(0.0, 1.0) >= mid: #***bug*** it should be <= (Fixed in 0.0.2)
                        selected_introns.append(i)
                        intron_prop_dict[i[0]] = float(np.random.choice(prop_list).strip())/100
                        num += 1
                elif i[-1] < 100:
                    if np.random.uniform(0.0, 1.0) >= bot: #***bug*** it should be <= (Fixed in 0.0.2)
                        selected_introns.append(i)
                        intron_prop_dict[i[0]] = float(np.random.choice(prop_list).strip())/100
                        num += 1
            if num >= num_cap:
                break
    #Make gene dictionary for introns
    intron_dict = OrderedDict()
    for i in selected_introns:
        try:
            intron_dict[i[0]].append(i)
        except:
            intron_dict[i[0]] = []
            intron_dict[i[0]].append(i)
    for gene in intron_dict:
        if strand_dict[gene] == '+':
            intron_dict[gene].sort(key=takeSeventh)
        elif strand_dict[gene] == '-':
            intron_dict[gene].sort(key=takeSeventh, reverse=True)
    return dict(intron_dict), dict(intron_prop_dict)

#Make transcripts
def make_transcripts(gene_list, intron_dict, transcript_dict, ref_file, cdna_file, strand_dict):
    #Load reference genome
    ref = Fasta(ref_file)
    #Load cDNA
    cdna = Fasta(cdna_file)
    fasta_dict = OrderedDict()
    cdna_size_dict = OrderedDict()
    intron_size_dict = OrderedDict()
    junc_dict = OrderedDict()
    #Calculate total length in bases
    total_len = 0
    for gene in gene_list:
        #Check if gene transcript is in cdna file
        if transcript_dict[gene] not in cdna:
            sys.exit("ERROR: Transcript %s sequence not found in cDNA FASTA file." % (transcript_dict[gene]))
        if gene in intron_dict:
            #Check if gene chromosome is in reference genome file
            if intron_dict[gene][0][2] not in ref:
                sys.exit("ERROR: Chromosome %s from gene %s not found in reference genome." % (intron_dict[gene][0][2], gene))
            #cdna_counts = fpkm_dict[gene]*(1-intron_dict[gene][0][-1]) #Use intronic transcript proportion from first intron entry in dictionary
            #intronic_counts = fpkm_dict[gene]-cdna_counts
            #Construct intronic transcript
            current_pos = 0
            intronic_fasta = ''
            if strand_dict[gene] == '+':
                for intron in intron_dict[gene]:
                    intronic_fasta = intronic_fasta + str(cdna[transcript_dict[gene]][current_pos:intron[5]+1]) + str(ref[intron[2]][intron[3]:intron[4]+1])
                    junc_dict[transcript_dict[gene] + '-ei' + str(intron[6])] = len(intronic_fasta)-intron[7]
                    junc_dict[transcript_dict[gene] + '-ie' + str(intron[6])] = len(intronic_fasta)
                    #junc_dict[transcript_dict[gene] + '-ee' + str(intron[6])] = intron[5]
                    current_pos = intron[5]+1
                intronic_fasta = intronic_fasta + str(cdna[transcript_dict[gene]][current_pos:])
                fasta_dict[transcript_dict[gene]+'-e'] = cdna[transcript_dict[gene]][:]
            elif strand_dict[gene] == '-':
                for intron in intron_dict[gene]:
                    intronic_fasta = intronic_fasta + str(rc(str(cdna[transcript_dict[gene]][:]))[current_pos:intron[5]+1]) + str(ref[intron[2]][intron[3]:intron[4]+1])
                    junc_dict[transcript_dict[gene] + '-ei' + str(intron[6])] = len(intronic_fasta)-intron[7]
                    junc_dict[transcript_dict[gene] + '-ie' + str(intron[6])] = len(intronic_fasta)
                    #junc_dict[transcript_dict[gene] + '-ee' + str(intron[6])] = intron[5]
                    current_pos = intron[5]+1
                intronic_fasta = intronic_fasta + str(rc(str(cdna[transcript_dict[gene]][:]))[current_pos:])
                fasta_dict[transcript_dict[gene]+'-e'] = rc(str(cdna[transcript_dict[gene]][:]))
            fasta_dict[transcript_dict[gene]+'-i'] = intronic_fasta
            intron_size_dict[transcript_dict[gene]] = len(intronic_fasta)
            total_len = total_len + len(cdna[transcript_dict[gene]][:])#*cdna_counts)
            total_len = total_len + len(intronic_fasta)#*intronic_counts)
            cdna_size_dict[transcript_dict[gene]] = len(cdna[transcript_dict[gene]][:])
        else:
            #cdna_counts = fpkm_dict[gene]
            if strand_dict[gene] == '+':
                fasta_dict[transcript_dict[gene]+'-e'] = cdna[transcript_dict[gene]][:]
            elif strand_dict[gene] == '-':
                fasta_dict[transcript_dict[gene]+'-e'] = rc(str(cdna[transcript_dict[gene]][:]))
            total_len = total_len + len(cdna[transcript_dict[gene]][:])#*cdna_counts)
            cdna_size_dict[transcript_dict[gene]] = len(cdna[transcript_dict[gene]][:])
    return dict(fasta_dict), total_len, dict(cdna_size_dict), dict(intron_size_dict), dict(junc_dict)

#Calculate and adjust FPKM model to fit number of aligned reads
def adjust_fpkm(gene_list, fpkm_dict, num_read, cdna_size_dict, transcript_dict, ran_read):
    #calculate total counts
    counts = 0
    for gene in gene_list:
        counts = counts + (fpkm_dict[gene]*cdna_size_dict[transcript_dict[gene]]*(num_read*(1-ran_read)/2))/1000000000
    #calculate ratio of actual counts to requested counts
    r = counts*2/(num_read*(1-ran_read))
    return r

#Coverage saturation checker (check that the minimum FPKM translate to at least 10 fragments for shortest cDNA length according to number of total reads)
def coverage_check(gene_list, fpkm_dict, lowest_FPKM, num_read, cdna_size_dict, transcript_dict, biotype_dict, ins_length, ins_length_stdev, ran_read, trans_gene_dict):
    #adjust FPKM
    correction_ratio = adjust_fpkm(gene_list, fpkm_dict, num_read, cdna_size_dict, transcript_dict, ran_read)
    new_fpkm_dict = OrderedDict()
    for gene in gene_list:
        new_fpkm_dict[gene] = fpkm_dict[gene]/correction_ratio
    corrected_lowest_FPKM = lowest_FPKM/correction_ratio
    #Get shortest cDNA length
    short_cdna = 1000000 #begin at 1000000 bp
    for transcript in cdna_size_dict:
        try: #In case transcript not found in GTF file
            if biotype_dict[trans_gene_dict[transcript]] == 'protein_coding':
                if short_cdna > cdna_size_dict[transcript] and cdna_size_dict[transcript] > (ins_length + (ins_length_stdev*3)): #ins_length + (ins_length_stdev*3) DWGSIM limit
                    short_cdna = cdna_size_dict[transcript]
        except:
            pass
    if int((corrected_lowest_FPKM*(short_cdna/1000))*((num_read*(1-ran_read)/2)/1000000)) < 10:
        #Calculate suggesting values
        suggest_num_read = math.ceil(((10/(corrected_lowest_FPKM*(short_cdna/1000)))*1000000)/(1-ran_read))*2
        suggest_min_FPKM = round((10/((num_read*(1-ran_read)/2)/1000000))/(short_cdna/1000)*correction_ratio,3)
        if suggest_num_read < 1000000000:
            sys.exit("ERROR: Insufficient coverage, please either increase number of reads to %s or increase the lower FPKM limit to %s." % (suggest_num_read, suggest_min_FPKM))
        else:
            sys.exit("ERROR: Insufficient coverage, please increase the lower FPKM limit to %s." % (suggest_min_FPKM))
    return dict(new_fpkm_dict)

#Check transcript lengths
def check_cdna_len(cdna_size_dict, ins_length, ins_length_stdev, fasta_dict):
    omit_gene_len_list = []
    check = 0
    for key in cdna_size_dict:
        if cdna_size_dict[key] <= (ins_length + 3*ins_length_stdev):
            check = 1
            omit_gene_len_list.append(key)
    if check == 1:
        print("INFO: Transcripts with length less than %s bp will not be simulated." % (str(round(ins_length + 3*ins_length_stdev))))
    for transcript in omit_gene_len_list:
        del fasta_dict[transcript + '-e']
        try:
            del fasta_dict[transcript + '-i']
        except:
            pass
    return fasta_dict, omit_gene_len_list

#Write transcriptome FASTA file
def write_fasta(fasta_dict, output_dir):
    file_path = os.path.join(output_dir, 'transcriptome.fa')
    out = open(file_path, 'w')
    for transcript in fasta_dict:
        out.write('{}\n{}\n'.format('>'+transcript, str(fasta_dict[transcript])))
    out.close()

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
    sim1_num_frags = int((transcript_frag_dict[last_rep_name]['@'+transcript_dict[gene_name]+'-e']*total_len)/cdna_size_dict[transcript_dict[gene_name]])*total_reps
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

#Clean fastq files
def cleaner(ref):
    check_call(['rm', ref + '.out.bwa.read1.fastq.gz'], stdout=DEVNULL, stderr=STDOUT)
    check_call(['rm', ref + '.out.bwa.read2.fastq.gz'], stdout=DEVNULL, stderr=STDOUT)
    check_call(['rm', ref + '.out.bwa1.pkl'], stdout=DEVNULL, stderr=STDOUT)
    check_call(['rm', ref + '.out.bwa2.pkl'], stdout=DEVNULL, stderr=STDOUT)
    check_call(['rm', ref], stdout=DEVNULL, stderr=STDOUT)
    check_call(['rm', ref + '.out.mutations.vcf'], stdout=DEVNULL, stderr=STDOUT)
    check_call(['rm', ref + '.out.mutations.txt'], stdout=DEVNULL, stderr=STDOUT)

#Clear memory and storage
def multi_clear(current_gene_list, fastq_dict1, fastq_dict2, output_dir, num_threads, rep_names_dict):
    fastq_dict1.clear()
    fastq_dict2.clear()
    in_fasta_list = []
    if len(current_gene_list) >= num_threads:
        for index in range(num_threads):
            in_fasta_list.append(os.path.join(output_dir, 'tmp' + str(index) + '.fa'))
    else:
        for index in range(len(current_gene_list)):
            in_fasta_list.append(os.path.join(output_dir, 'tmp' + str(index) + '.fa'))
    jobs = []
    for ref in in_fasta_list:
        p = multiprocessing.Process(target=cleaner, args=(ref,))
        jobs.append(p)
        p.start()
    for rep in rep_names_dict:
        try:
            check_call(['rm', os.path.join(output_dir, rep_names_dict[rep]) + '1.pkl'], stdout=DEVNULL, stderr=STDOUT)
            check_call(['rm', os.path.join(output_dir, rep_names_dict[rep]) + '2.pkl'], stdout=DEVNULL, stderr=STDOUT)
        except:
            pass
    for job in jobs:
        job.join()

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
                        fastq_dict1[line.rsplit('_', 9)[0]].append([line, next(f1), next(f1), next(f1)])
                    except:
                        fastq_dict1[line.rsplit('_', 9)[0]] = []
                        fastq_dict1[line.rsplit('_', 9)[0]].append([line, next(f1), next(f1), next(f1)])
                #Elif read1 is reverse
                elif line.rsplit('_', 7)[1] == '1':
                    #Place in dict2 and change read number
                    try:
                        fastq_dict2[line.rsplit('_', 9)[0]].append([line.rsplit('/', 1)[0] + '/2\n', next(f1), next(f1), next(f1)])
                    except:
                        fastq_dict2[line.rsplit('_', 9)[0]] = []
                        fastq_dict2[line.rsplit('_', 9)[0]].append([line.rsplit('/', 1)[0] + '/2\n', next(f1), next(f1), next(f1)])
            elif strand_dict[trans_gene_dict[line.rsplit('_', 9)[0][1:-2]]] == '-':
                #If read1 is forward
                if line.rsplit('_', 7)[1] == '0':
                    #Place in dict2 and change read number
                    try:
                        fastq_dict2[line.rsplit('_', 9)[0]].append([line.rsplit('/', 1)[0] + '/2\n', next(f1), next(f1), next(f1)])
                    except:
                        fastq_dict2[line.rsplit('_', 9)[0]] = []
                        fastq_dict2[line.rsplit('_', 9)[0]].append([line.rsplit('/', 1)[0] + '/2\n', next(f1), next(f1), next(f1)])
                #Elif read1 is reverse
                elif line.rsplit('_', 7)[1] == '1':
                    #Place in dict1
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
                        fastq_dict2[line.rsplit('_', 9)[0]].append([line, next(f2), next(f2), next(f2)])
                    except:
                        fastq_dict2[line.rsplit('_', 9)[0]] = []
                        fastq_dict2[line.rsplit('_', 9)[0]].append([line, next(f2), next(f2), next(f2)])
                #Elif read1 is reverse
                elif line.rsplit('_', 7)[1] == '1':
                    #Place in dict1 and change read number
                    try:
                        fastq_dict1[line.rsplit('_', 9)[0]].append([line.rsplit('/', 1)[0] + '/1\n', next(f2), next(f2), next(f2)])
                    except:
                        fastq_dict1[line.rsplit('_', 9)[0]] = []
                        fastq_dict1[line.rsplit('_', 9)[0]].append([line.rsplit('/', 1)[0] + '/1\n', next(f2), next(f2), next(f2)])
            elif strand_dict[trans_gene_dict[line.rsplit('_', 9)[0][1:-2]]] == '-':
                #If read1 is forward
                if line.rsplit('_', 7)[1] == '0':
                    #Place in dict1 and change read number
                    try:
                        fastq_dict1[line.rsplit('_', 9)[0]].append([line.rsplit('/', 1)[0] + '/1\n', next(f2), next(f2), next(f2)])
                    except:
                        fastq_dict1[line.rsplit('_', 9)[0]] = []
                        fastq_dict1[line.rsplit('_', 9)[0]].append([line.rsplit('/', 1)[0] + '/1\n', next(f2), next(f2), next(f2)])
                #Elif read1 is reverse
                elif line.rsplit('_', 7)[1] == '1':
                    #Place in dict2
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

#Multiprocessing extract_fastq (Pickle) stranded
def multi_thread_extract_fastq_pickle(current_gene_list, output_dir, trans_gene_dict, strand_dict, num_threads):    
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

#Make output files
def make_out_files(output_dir, transcript_frag_dict, intron_reps, control_reps):
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
            #try:
            chunks_dict[rep][transcript_id] = [latest_index[transcript_id], min(transcript_frag_dict[rep][transcript_id] + latest_index[transcript_id], len(fastq_dict1[transcript_id]))]
            latest_index[transcript_id] = min(transcript_frag_dict[rep][transcript_id] + latest_index[transcript_id], len(fastq_dict1[transcript_id]))
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

#Create output info file
def output_info(intron_list, intron_dict, real_frag_dict, read_junc_dict, output_dir, strand_dict, fpkm_dict, intron_prop_dict, rep_names_dict, omit_gene_len_list):
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
        if intron[1] not in omit_gene_len_list:
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

def main():
    if len(argv) != 6:
        sys.exit("Usage: irsim-0.0.1-pre-alpha.py ref_genome.fa cDNA.fa annotation.gtf FPKM_model.tsv config.ini")
    now = datetime.datetime.now()
    print("RNA IR Simulation started...")
    print(now)
    ref_file = argv[1]
    cdna_file = argv[2]
    gtf_file = argv[3]
    model_file = argv[4]
    config_file = argv[5]
    start_time = time.time()
    #Setup config
    config = configparser.ConfigParser()
    config.read(config_file)
    #Check config.ini file
    config_checker(config)
    output_dir = config['Paths']['Output directory']
    dwgsim_dir = config['Paths']['DWGSIM directory']
    intron_reps = config['Sequencing-details']['Number of replicates for dataset with intron retention']
    control_reps = config['Sequencing-details']['Number of replicates for dataset without intron retention']
    total_reps = int(intron_reps) + int(control_reps)
    max_rep_var = config['Sequencing-details']['Maximum percent variance between replicates']
    num_read = config['Sequencing-details']['Number of reads per replicate']
    #paired = config['Sequencing-details']['Paired-end sequencing']
    stranded = config['Sequencing-details']['Strand-specific']
    read_length = config['Sequencing-details']['Read length']
    ins_length = config['Sequencing-details']['Average outer distance between read pairs (insert length)']
    ins_length_stdev = config['Sequencing-details']['Standard deviation of outer distance between read pairs']
    first_read_error = config['Sequencing-details']['Per base error rate of the first read']
    second_read_error = config['Sequencing-details']['Per base error rate of the second read']
    mut_rate = config['Sequencing-details']['Rate of mutation']
    indel_rate = config['Sequencing-details']['Fraction of mutations that are indels']
    indel_ext = config['Sequencing-details']['Probability an indel is extended']
    indel_len = config['Sequencing-details']['Minimum length indel']
    ran_read = config['Sequencing-details']['Probability of a random read']
    mock_ran_read = 0
    perc_intron = config['Intron-retention']['Percentage of total introns which are retained']
    intron_len1 = config['Intron-retention']['Probability of retaining an intron with length more than 1000 bases']
    intron_len2 = config['Intron-retention']['Probability of retaining an intron with length between 100-1000 bases']
    intron_len3 = config['Intron-retention']['Probability of retaining an intron with length less than 100 bases']
    intron_proportion_list = config['Intron-retention']['Proportions of transcripts with intron retention']
    lower_FPKM = config['Gene-expression']['Lower limit of FPKM in gene expression model']
    upper_FPKM = config['Gene-expression']['Upper limit of FPKM in gene expression model']
    seed = config['Seed']['Random seed']
    num_threads = config['Threads']['Number of threads']
    #Assign random seed (For same intron simulation)
    assign_seed(seed)
    #Check output_dir exist
    out_dir_exist(output_dir)
    #Gene and number dictionary generation from FPKM_model
    fpkm_dict, lowest_FPKM, highest_FPKM = ge_extract(model_file, float(lower_FPKM), float(upper_FPKM))
    #Obtain intron information of selected genes
    gene_list, intron_list, transcript_dict, strand_dict, biotype_dict, trans_gene_dict, transcript_coord_dict = get_intron_info(gtf_file, fpkm_dict)
    #Select introns to retain
    intron_dict, intron_prop_dict = select_introns(intron_list, float(perc_intron)/100, float(intron_len1), float(intron_len2), float(intron_len3), intron_proportion_list.split(','), strand_dict)
    #Make intronic transcripts
    fasta_dict, total_len, cdna_size_dict, intron_size_dict, junc_dict = make_transcripts(gene_list, intron_dict, transcript_dict, ref_file, cdna_file, strand_dict)
    #Check transcript lengths
    fasta_dict, omit_gene_len_list = check_cdna_len(cdna_size_dict, float(ins_length), float(ins_length_stdev), fasta_dict)
    #Coverage saturation checker
    fpkm_dict = coverage_check(gene_list, fpkm_dict, lowest_FPKM, int(num_read), cdna_size_dict, transcript_dict, biotype_dict, int(ins_length), int(ins_length_stdev), float(ran_read), trans_gene_dict)
    #Write transcriptome FASTA file 1
    write_fasta(fasta_dict, output_dir)
    #Decipher amount of simulated reads to suit 20% of transcripts FPKM
    sim1_num_frags, percentile_FPKM = sim_amount(gene_list, fpkm_dict, int(num_read), total_len, 20, fasta_dict, total_reps)
    print("--- %s seconds --- Simulation part 1.0/5.0 configuration" % (time.time() - start_time))
    #Split reference fasta to multiple fasta files according to number of threads
    split_fasta(output_dir, 'transcriptome.fa', int(num_threads), len(fasta_dict))
    #Run multi dwgsim
    multi_subprocess_dwgsim(gene_list, dwgsim_dir, output_dir, int(num_threads), first_read_error, second_read_error, ins_length, ins_length_stdev, sim1_num_frags, read_length, mut_rate, indel_rate, indel_ext, indel_len, mock_ran_read, seed)
    print("--- %s seconds --- Simulation part 1.0/5.0 reads generated" % (time.time() - start_time))
    #Create list for random reads
    rand_read_list = []
    #Multiprocessing extract fastq
    fastq_dict1, fastq_dict2 = multi_thread_extract_fastq_pickle(gene_list, output_dir, trans_gene_dict, strand_dict, int(num_threads))
    print("--- %s seconds --- Simulation part 1.0/5.0 reads extracted" % (time.time() - start_time))
    #First and only make transcript_frag_dict for each rep
    transcript_frag_dict = OrderedDict()
    real_frag_dict = OrderedDict()
    for i in range(int(intron_reps)):
        transcript_frag_dict['i-' + str(len(transcript_frag_dict))], real_frag_dict['i-' + str(len(real_frag_dict))] = make_fragment_dict(fastq_dict1, intron_dict, intron_prop_dict, fpkm_dict, trans_gene_dict, int(num_read), cdna_size_dict, float(max_rep_var), 'intron', float(ran_read), intron_size_dict)
    #spacer
    for i in range(int(control_reps)):
        transcript_frag_dict['c-' + str(len(transcript_frag_dict))], real_frag_dict['c-' + str(len(real_frag_dict))] = make_fragment_dict(fastq_dict1, intron_dict, intron_prop_dict, fpkm_dict, trans_gene_dict, int(num_read), cdna_size_dict, float(max_rep_var), 'control', float(ran_read), intron_size_dict)
    #spacer
    #Make out files
    rep_names_dict = make_out_files(output_dir, transcript_frag_dict, int(intron_reps), int(control_reps))
    #Make read index chunks for each transcript
    chunks_dict = read_index_chunks(transcript_frag_dict, fastq_dict1)
    #Prepare dicts for junc read counting
    read_junc_dict, total_intron_dict = junc_count_prep(intron_list, junc_dict, rep_names_dict)
    #Create fastq lists and write fastq
    multi_fastq_list_write(num_threads, transcript_frag_dict, fastq_dict1, fastq_dict2, total_reps, output_dir, rep_names_dict, chunks_dict, junc_dict, read_junc_dict, total_intron_dict, read_length)
    print("--- %s seconds --- Simulation part 1.0/5.0 completed" % (time.time() - start_time))
    #Clear memory and delete files
    multi_clear(gene_list, fastq_dict1, fastq_dict2, output_dir, int(num_threads), rep_names_dict)
    #Multiple rounds with increasing percentile
    for percentile in [40, 60, 80, 90, 100]:
        new_len, new_gene_list = update_fasta_dict(fasta_dict, transcript_frag_dict, trans_gene_dict)
        write_fasta(fasta_dict, output_dir)
        sim1_num_frags, percentile_FPKM = sim_amount2(new_gene_list, fpkm_dict, new_len, percentile, total_reps, transcript_frag_dict, transcript_dict, cdna_size_dict)
        print("--- %s seconds --- Simulation part %s/5.0 configuration" % (time.time() - start_time, round(percentile/20,1)))
        split_fasta(output_dir, 'transcriptome.fa', int(num_threads), len(fasta_dict))
        multi_subprocess_dwgsim(new_gene_list, dwgsim_dir, output_dir, int(num_threads), first_read_error, second_read_error, ins_length, ins_length_stdev, sim1_num_frags, read_length, mut_rate, indel_rate, indel_ext, indel_len, mock_ran_read, seed)
        print("--- %s seconds --- Simulation part %s/5.0 reads generated" % (time.time() - start_time, round(percentile/20,1)))
        fastq_dict1, fastq_dict2 = multi_thread_extract_fastq_pickle(new_gene_list, output_dir, trans_gene_dict, strand_dict, int(num_threads))
        print("--- %s seconds --- Simulation part %s/5.0 reads extracted" % (time.time() - start_time, round(percentile/20,1)))
        chunks_dict = read_index_chunks(transcript_frag_dict, fastq_dict1)
        multi_fastq_list_write(num_threads, transcript_frag_dict, fastq_dict1, fastq_dict2, total_reps, output_dir, rep_names_dict, chunks_dict, junc_dict, read_junc_dict, total_intron_dict, read_length)
        print("--- %s seconds --- Simulation part %s/5.0 completed" % (time.time() - start_time, round(percentile/20,1)))
        multi_clear(new_gene_list, fastq_dict1, fastq_dict2, output_dir, int(num_threads), rep_names_dict)
    #spacer
    percentile = 100
    new_len, new_gene_list = update_fasta_dict(fasta_dict, transcript_frag_dict, trans_gene_dict)
    while len(fasta_dict) > 0:
        write_fasta(fasta_dict, output_dir)
        sim1_num_frags, percentile_FPKM = sim_amount2(new_gene_list, fpkm_dict, new_len, percentile, total_reps, transcript_frag_dict, transcript_dict, cdna_size_dict)
        print("--- %s seconds --- Simulation part %s/5.0 configuration (Extra)" % (time.time() - start_time, round(percentile/20,1)))
        split_fasta(output_dir, 'transcriptome.fa', int(num_threads), len(fasta_dict))
        multi_subprocess_dwgsim(new_gene_list, dwgsim_dir, output_dir, int(num_threads), first_read_error, second_read_error, ins_length, ins_length_stdev, sim1_num_frags, read_length, mut_rate, indel_rate, indel_ext, indel_len, mock_ran_read, seed)
        print("--- %s seconds --- Simulation part %s/5.0 reads generated (Extra)" % (time.time() - start_time, round(percentile/20,1)))
        fastq_dict1, fastq_dict2 = multi_thread_extract_fastq_pickle(new_gene_list, output_dir, trans_gene_dict, strand_dict, int(num_threads))
        print("--- %s seconds --- Simulation part %s/5.0 reads extracted (Extra)" % (time.time() - start_time, round(percentile/20,1)))
        chunks_dict = read_index_chunks(transcript_frag_dict, fastq_dict1)
        multi_fastq_list_write(num_threads, transcript_frag_dict, fastq_dict1, fastq_dict2, total_reps, output_dir, rep_names_dict, chunks_dict, junc_dict, read_junc_dict, total_intron_dict, read_length)
        print("--- %s seconds --- Simulation part %s/5.0 completed (Extra)" % (time.time() - start_time, round(percentile/20,1)))
        multi_clear(new_gene_list, fastq_dict1, fastq_dict2, output_dir, int(num_threads), rep_names_dict)
    #spacer
    #Simulate random reads
    #Count existing reads of each dataset
    read_count_dict = read_count(output_dir, rep_names_dict)
    num_ran_reads = random_read_counter(read_count_dict, int(num_read), float(ran_read))
    #Create mock input fasta
    mock_fasta(output_dir, int(num_threads))
    #Simulate random reads
    multi_subprocess_dwgsim(gene_list, dwgsim_dir, output_dir, int(num_threads), first_read_error, second_read_error, ins_length, ins_length_stdev, num_ran_reads, read_length, mut_rate, indel_rate, indel_ext, indel_len, 1, seed)
    print("--- %s seconds --- Simulated random reads" % (time.time() - start_time))
    ran_fastq_dict1, ran_fastq_dict2 = multi_thread_extract_random_fastq_pickle(output_dir, trans_gene_dict, strand_dict, int(num_threads))
    random_reads_write(read_count_dict, ran_fastq_dict1, ran_fastq_dict2, output_dir, rep_names_dict, int(num_read))
    print("--- %s seconds --- Finished writing random reads." % (time.time() - start_time))
    #Clear memory and delete files
    multi_clear(gene_list, fastq_dict1, fastq_dict2, output_dir, int(num_threads), rep_names_dict)
    check_call(['rm', os.path.join(output_dir, 'transcriptome.fa')], stdout=DEVNULL, stderr=STDOUT)
    print("--- %s seconds --- Prepared report information." % (time.time() - start_time))
    #Writing output report file
    output_info(intron_list, intron_dict, real_frag_dict, read_junc_dict, output_dir, strand_dict, fpkm_dict, intron_prop_dict, rep_names_dict, omit_gene_len_list)
    #Compress fastq to gz format
    compress_multi(rep_names_dict, output_dir, int(num_threads), total_reps)
    print("--- %s seconds --- RNA IR Simulation finished." % (time.time() - start_time))
    now = datetime.datetime.now()
    print(now)

if __name__ == '__main__':
    main()


