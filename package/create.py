"""
Creator.py

This module simulates a transcriptome and randomly creates intron retained transcripts.

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

import sys
import math
import os.path
import numpy as np
from pyfaidx import Fasta
from collections import OrderedDict

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

#Parse GTF attribute column into dictionary, also does checking
def gtf_attribute_parser_gene_id(gtf_line, gtf_line_next, index, index_next):
    attr_dict = OrderedDict()
    attr_dict_next = OrderedDict()
    attribute_list = gtf_line.split('\t')[8].split(';')
    attribute_list_next = gtf_line_next.split('\t')[8].split(';')
    attribute_list = filter(None, attribute_list)
    attribute_list_next = filter(None, attribute_list_next)
    for item in attribute_list:
        attr_dict[item.strip().split(' ')[0]] = item.strip().split(' ')[1].strip('\"\'')
    if 'gene_id' not in attr_dict:
        sys.exit("ERROR: gene_id not found in GTF attributes column in line %s." % str(index+1))
    for item in attribute_list_next:
        attr_dict_next[item.strip().split(' ')[0]] = item.strip().split(' ')[1].strip('\"\'')
    if 'gene_id' not in attr_dict_next:
        sys.exit("ERROR: gene_id not found in GTF attributes column in line %s." % str(index_next+1))
    return dict(attr_dict), dict(attr_dict_next)

#Parse GTF attribute column into dictionary, also does checking
def gtf_attribute_parser_exon(gtf_line, index):
    attr_dict = OrderedDict()
    attribute_list = gtf_line.split('\t')[8].split(';')
    attribute_list = filter(None, attribute_list)
    for item in attribute_list:
        attr_dict[item.strip().split(' ')[0]] = item.strip().split(' ')[1].strip('\"\'')
    if 'gene_id' not in attr_dict:
        sys.exit("ERROR: gene_id not found in GTF attributes column in line %s." % str(index+1))
    if 'transcript_id' not in attr_dict:
        sys.exit("ERROR: transcript_id not found in GTF attributes column in line %s." % str(index+1))
    if 'gene_biotype' not in attr_dict:
        if 'gene_type' not in attr_dict:
            sys.exit("ERROR: gene_biotype/gene_type not found in GTF attributes column in line %s." % str(index+1))
    else:
        attr_dict['gene_type'] = attr_dict['gene_biotype']
    return dict(attr_dict)

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
            attr_dict, attr_dict_next = gtf_attribute_parser_gene_id(gtf[index], gtf[index+1], index, index+1)
            #Check next line for same gene_id and make exon list(Ensure GTF uses double quotes and correct columns) Col: gene_id,transcript_id,chr,start,end,strand
            if attr_dict['gene_id'] == attr_dict_next['gene_id']:
                if gtf[index].split('\t')[2] == "exon":
                    attr_dict = gtf_attribute_parser_exon(gtf[index], index)
                    exon_tmp_list.append([attr_dict['gene_id'], attr_dict['transcript_id'], str(gtf[index].split('\t')[0]), int(gtf[index].split('\t')[3]), int(gtf[index].split('\t')[4]), str(gtf[index].split('\t')[6])])
                    biotype_dict[attr_dict['gene_id']] = attr_dict['gene_type']
            else:
                if gtf[index].split('\t')[2] == "exon":
                    attr_dict = gtf_attribute_parser_exon(gtf[index], index)
                    exon_tmp_list.append([attr_dict['gene_id'], attr_dict['transcript_id'], str(gtf[index].split('\t')[0]), int(gtf[index].split('\t')[3]), int(gtf[index].split('\t')[4]), str(gtf[index].split('\t')[6])])
                    biotype_dict[attr_dict['gene_id']] = attr_dict['gene_type']
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
                                intron_list.append([exon_tmp_list[index][0], exon_tmp_list[index][1], exon_tmp_list[index][2], int(exon_tmp_list[index][4]+1), int(exon_tmp_list[index+1][3]-1), int(exon_tmp_list[index][4]+1-tsc-intron_size), index+1, int(exon_tmp_list[index+1][3]-1-(exon_tmp_list[index][4]+1)+1), ])
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
                    if np.random.uniform(0.0, 1.0) < top:
                        selected_introns.append(i)
                        intron_prop_dict[i[0]] = float(np.random.choice(prop_list).strip())/100
                        num += 1
                elif 100 <= i[-1] <= 1000:
                    if np.random.uniform(0.0, 1.0) < mid:
                        selected_introns.append(i)
                        intron_prop_dict[i[0]] = float(np.random.choice(prop_list).strip())/100
                        num += 1
                elif i[-1] < 100:
                    if np.random.uniform(0.0, 1.0) < bot:
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

#Generate on intron retention proportions for control sample
def select_ctrl_intron_prop(intron_dict, prop_list):
    intron_prop_dict_ctrl = OrderedDict()
    for intron in intron_dict:
        intron_prop_dict_ctrl[intron] = float(np.random.choice(prop_list).strip())/100
    return dict(intron_prop_dict_ctrl)

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
            sys.exit("ERROR: Insufficient coverage, please either increase minimum number of reads to %s or increase the lower FPKM limit to %s." % (suggest_num_read, suggest_min_FPKM))
        else:
            sys.exit("ERROR: Insufficient coverage, please increase the lower FPKM limit to %s." % (suggest_min_FPKM))
        #print("ERROR: Insufficient coverage, please either increase number of reads to %s or increase the lower FPKM limit to %s." % (suggest_num_read, suggest_min_FPKM))
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
