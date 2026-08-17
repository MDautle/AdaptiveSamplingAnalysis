# Adaptive Sampling Analysis

This repository is associated with the following manuscript:  Dautle, M., Warren, M., Roman, E., Gould, A. (2026) Evaluating Adaptive Sampling Performance for Low-Diversity Whole-Genome Sequencing Using a Binary Vertebrate-Bacteria Symbiosis. (Submitted)

## Overview 
This repository contains the scripts and source data required to reproduce the figures in the associated paper. The repository includes the following files: 
| File | Description |
|------|-------------|
| [AdaptiveSamplingAssembled_v06-15-2026.ipynb](./AdaptiveSamplingAssembled_v06-15-2026.ipynb) | The main set of scripts used to produce source data from raw data and associated visualizations. Minimal changes to redact file paths have been made, but it is otherwise unaltered. All packages are listed in cell [1]. Some packages have syntax changes between versions and can generate errors when not set to a specific version. Those packages are detailed with the version number in cell [2].  | 
| [AdaptiveSampling_SourceData.xlsb](./AdaptiveSampling_SourceData.xlsb) | The source data for all figure panels. Each tab is named for the figure panel(s) the data are associated with. These source data were generated from this script; however, column names and formatting have been updated for clarity and to reduce the file size. The original, unaltered source data files are available upon request. | 
| [plot_bedMethyl_HeatmapTrack_withGATC.py](./plot_bedMethyl_HeatmapTrack_withGATC.py) | A helper script for the Jupyter Notebook. It is only required for Figure 3C. |

## Reproducing the Figures
If you wish to reproduce the figures, a few changes will be required. Those changes are as follows: 
+ The file paths must be updated to match the file locations on your device
+ Column names may have been changed on the source data file during compilation. 
+ The source data files were compiled into one Excel sheet to reduce file size requirements for long-term storage. However, the notebook is currently designed for each tab to be read in as and independent csv file (with the exception of Figure 1A, each column is its own file). You will either need to change how the file is read in (i.e. read in the [source data file](./AdaptiveSampling_SourceData.xlsb) as a dictionary of pandas dataframes with pd.read_excel), or save each tab as an independent file.

As a reminder, not all cells will need to be run when generating the figures from the source data. The script contains other processing steps which are not necessary when starting with the source data files. 


## Minimally Processed Raw Data 
The ONT long-read sequencing data in FASTQ format can be found on NCBI under BioProject ID PRJNA1504951. It will be released upon manuscript publication. If you would like early access, you may make a request by emailing Ms.Madison Dautle (madison.dautle@temple.edu). Please include the reason for your request in your email. As FASTQ format does not retain methylation modification tags, the files are available in modBam format or pod5 format upon reasonable request. 

## Contact Us
The best way to reach us is by email, rather than opening an issue on the repository. Questions may be directed to Dr.Alison Gould (alison.gould@temple.edu) or Ms.Madison Dautle (madison.dautle@temple.edu). 
