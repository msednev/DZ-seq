# High-throughput activity analysis of RNA-cleaving deoxyribozymes by DZ-seq

Source code for publication **High-throughput activity analysis of RNA-cleaving deoxyribozymes by DZ-seq** 
by Maksim V. Sednev, Anam Liaqat and Claudia Höbartner

##### Project structure

```
DZ-seq
├── input                         <- Input files
│   └── raw_reads                 <- Raw sequencing reads
├── notebooks                     <- Jupyter notebooks for analysis of generated data and plotting figures
├── output                        <- Data generated during processing steps
│   ├── counts                    <- Counts of unique deoxyribozyme+RNA combinations for each library
│   ├── datasets                  <- Datasets on counts and cleavage activity
│   └── processed_reads           <- Sequencing reads trimmed by cutadapt
│       ├── core                  <- Sequencing reads trimmed down to deoxyribozymes cores
│       │   └── logs              <- cutadapt logs
│       └── substrate             <- Sequencing reads trimmed down to RNA substrate residues
│           └── logs              <- cutadapt logs
└── scripts                       <- R scripts for data processing
```

#### Requirements
The DZ-seq pipeline was tested on Ubuntu 18.04 with the following packages:
- Python 3.8.8
- cutadapt 3.4
- numpy 1.18.1
- pandas 1.0.1
- matplotlib 3.3.4
- seaborn 0.11.2
- R 3.6.3
- data.table 1.14.0
- stringi 1.7.3