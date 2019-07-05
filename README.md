# IRSim: Simulation of Intron Retention in Coding RNA 

## (WORK IN PROGRESS)

Latest version 0.0.1-pre-alpha

This software is created to simulate RNA sequencing datasets (Illumina NGS) with pseudo-random intron retention events in coding RNA transcripts. The simulated datasets can be used to evaluate and benchmark intron-retention detection softwares or workflows. 

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
chmod +x irsim-0.0.1-pre-alpha.py
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
DWGSIM directory = /path/to/DWGSIM_dir
```
* Add path to Ouput directory
```
Output directory = /path/to/output_directory
```
* Edit all other parameters to your choice
```
Number of replicates for dataset with intron retention = 1 
Number of replicates for dataset without intron retention = 1 
.
.
.
Number of threads = 10
```
Please note that this version only allows maximum 1 replicate of each dataset due to high memory consumption.

### 3. Run IRSim
```
/path/to/irsim-0.0.1-pre-alpha.py ref_genome.fa cDNA.fa annotation.gtf FPKM_model.tsv config.ini
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
