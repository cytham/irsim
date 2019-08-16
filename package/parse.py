"""
Parse.py

This module checks for input files and their format.

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
import os.path
from collections import OrderedDict

#Check sections and keys of config.ini file
def config_checker(config):
    for key in ['DWGSIM directory']:
        if config.has_option('Paths', key) == False:
            print('Error in config.ini: "%s" key do not exists' % (key))
            sys.exit()
    for key in ['Number of replicates for sample A (Experimental sample)', 'Number of replicates for sample B (Control sample)', 'Simulate intron retention in sample B (yes/no)', 'Maximum percent variance of read counts between replicates (%)', 'Minimum number of reads per sample/replicate', 'Paired-end sequencing', 'Strand-specific (yes/no)', 'Read length in bases', 'Average outer distance between read pairs (aka insert length)', 'Standard deviation of outer distance between read pairs', "Per base error rate of first read (Range from 5' to 3' ends)", "Per base error rate of second read (Range from 5' to 3' ends)", 'Rate of mutation', 'Fraction of mutations that are indels', 'Probability that an indel is extended', 'Minimum length of indel', 'Probability of a random read']:
        if config.has_option('Sequencing details', key) == False:
            print('Error in config.ini: "%s" key do not exists' % (key))
            sys.exit()
    for key in ['Percentage of total introns which are retained (%)', 'Probability of retaining an intron with length more than 1000 bases', 'Probability of retaining an intron with length between 100-1000 bases', 'Probability of retaining an intron with length less than 100 bases', 'Set of possible Percent Intron Retention (PIR) for each intron (%)']:
        if config.has_option('Intron retention', key) == False:
            print('Error in config.ini: "%s" key do not exists' % (key))
            sys.exit()
    for key in ['Lower limit of FPKM in the provided gene expression model', 'Upper limit of FPKM in the provided gene expression model']:
        if config.has_option('Gene expression', key) == False:
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
    for key in ['Verbose (yes/no)']:
        if config.has_option('Verbosity', key) == False:
            print('Error in config.ini: "%s" key do not exists' % (key))
            sys.exit()

    #for key in ['Memory consumption level (1 ~ <50Gb; 2 ~ <100Gb; 3 ~ <150Gb)']:
        #if config.has_option('Memory consumption', key) == False:
            #print('Error in config.ini: "%s" key do not exists' % (key))
            #sys.exit()

#Check output dir exist
def check_dir_file_exist(ref_file, cdna_file, gtf_file, model_file, config_file, output_dir):
    if not os.path.exists(output_dir):
        sys.exit("ERROR: Output directory '%s' is not found." % (output_dir))
    if not os.path.isfile(ref_file):
        sys.exit("ERROR: The file '%s' is not found." % (ref_file))
    if not os.path.isfile(cdna_file):
        sys.exit("ERROR: The file '%s' is not found." % (cdna_file))
    if not os.path.isfile(gtf_file):
        sys.exit("ERROR: The file '%s' is not found." % (gtf_file))
    if not os.path.isfile(model_file):
        sys.exit("ERROR: The file '%s' is not found." % (model_file))
    if not os.path.isfile(config_file):
        sys.exit("ERROR: The file '%s' is not found." % (config_file))

#Configure replicates
def rep_config(total_intron_reps, total_control_reps):
    #rep_list = []
    total_rep_names_dict = OrderedDict()
    for index in range(int(total_intron_reps)):
        #rep_list.append("i-" + str(index))
        total_rep_names_dict["A-" + str(index)] = 'sample_A_rep' + str(index+1)
    for index in range(int(total_control_reps)):
        #rep_list.append("c-" + str(index))
        total_rep_names_dict["B-" + str(index)] = 'sample_B_rep' + str(index+1)
    return dict(total_rep_names_dict)
