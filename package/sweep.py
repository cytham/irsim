"""
Sweep.py

This module clears memory and removes unused files.

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
import multiprocessing
from subprocess import DEVNULL, STDOUT, check_call

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
