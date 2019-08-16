# IRSim: Simulation of Intron Retention in Coding RNA 

Latest version 0.0.2

This software is created to simulate RNA sequencing datasets (Illumina NGS) with pseudo-random intron retention events in coding RNA transcripts. The simulated datasets can be used to evaluate and benchmark intron-retention detection softwares or workflows. IRSim is a software package written in Python3 and employs DWGSIM for NGS read simulation.

The percent intron retention (PIR) of each retained intron is calculated by \[100 x mean retention reads divided by the sum of retention reads and spliced intron reads\] as seen in [Braunschweig et al., 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4216919/).

## Installation

### Requirements
1. Python 3 (3.6.7 or above)
2. numpy (1.16.2  or above) [Link](https://scipy.org/install.html)
3. pyfaidx (0.5.5.2  or above) [Link](https://pypi.org/project/pyfaidx/)
4. DWGSIM (0.1.11  or above) [Link](https://github.com/nh13/DWGSIM)

### Clone git repository:
```
git clone https://github.com/cytham/irsim.git 
cd irsim
chmod +x irsim
```

## Quick start
### 1. Required input files
Make sure you have the following input files ready:
1. A reference genome file in FASTA format
2. cDNA FASTA file
3. Annotation GTF file
4. A tab-delimited FPKM model file (Column 1 - Gene id \t Column 2 - FPKM values)

### 2. Edit config.ini file
* Add path to DWGSIM directory
```
DWGSIM directory = /path/to/DWGSIM_directory 
```
* Edit all other parameters to your choice
```
Number of threads = 10 
Number of replicates for sample A (Experimental sample) = 1
Number of replicates for sample B (Control sample) = 1
.
.
.
```

### 3. Run IRSim
```
/path/to/irsim ref_genome.fa cDNA.fa annotation.gtf FPKM_model.tsv config.ini output_directory
```

### 3. Output files
* Gzipped FASTQ paired-read files for each sample/replicate.
* A report file showing the Percent Intron Retention (PIR) for each intron in each sample/replicate.

## Versioning
See [Releases](https://github.com/cytham/irsim/releases)

## Citation

## Author

* **Tham Cheng Yong** - [cytham](https://github.com/cytham)

## License

This project is licensed under GNU General Public License - see [LICENSE](https://github.com/cytham/irsim/blob/master/LICENSE) for details.
