# IRSim: Simulation of Intron Retention in Coding RNA 

## (WORK IN PROGRESS)

Latest version 0.0.2

This software is created to simulate RNA sequencing datasets (Illumina NGS) with pseudo-random intron retention events in coding RNA transcripts. The simulated datasets can be used to evaluate and benchmark intron-retention detection softwares or workflows. 

## Installation

### Requirements
1. Python 3 (3.6.7 or above)
2. numpy (1.16.2  or above) [Link](https://scipy.org/install.html)
3. pyfaidx (0.5.5.2  or above) [Link](https://pypi.org/project/pyfaidx/)
4. DWGSIM (0.1.11  or above) [Link](https://github.com/nh13/DWGSIM)
5. A reference genome file in FASTA format
6. cDNA FASTA file
7. Annotation GTF file
8. An tab-delimited FPKM model file (Column 1 - Gene id \t Column 2 - FPKM values)

### Clone git repository:
```
git clone https://github.com/cytham/irsim.git 
cd irsim
chmod +x irsim
```

## Quick start
### 1. Edit config.ini file
* Add path to DWGSIM directory
```
DWGSIM directory = /path/to/DWGSIM_dir
```
* Edit all other parameters to your choice
```
Number of replicates for sample WITH intron retention (Experimental sample) = 2 
Number of replicates for sample WITHOUT intron retention (Control sample) = 2
.
.
.
Number of threads = 30
```

### 2. Run IRSim
```
/path/to/irsim ref_genome.fa cDNA.fa annotation.gtf FPKM_model.tsv config.ini output_directory
```

### 3. Output files
* Gzipped FASTQ paired-read files for each sample/replicate.
* A report file showing the Percent Intron Retention (PIR) for each intron in each sample/replicate.

## Versioning

## Citation

## Author

* **Tham Cheng Yong** - [cytham](https://github.com/cytham)

## License

This project is licensed under GNU General Public License - see [LICENSE](https://github.com/cytham/irsim/blob/master/LICENSE) for details.
