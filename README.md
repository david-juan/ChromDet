# ChromDet
A method to classify different cell types from their epigenomic profiles and detect genomic regions determining this classification.
This repository contains a modified version of S3det (Rausell et al. PNAS 2010) adapted to perform epigenomic analyses of sample segregation into different cell types. It contains several scripts aimed to facilitate this analysis by providing an interface with S3det and retrieving understable output files containing the Chromatin Sample Space and the Chromatin Determining Regions (CDRs).

## Installation

ChromDet has been tested on UNIX/Linux and MacOS systems.

### Directories

-  S3det distribution slightly modified to allow a higher number of variables and provide a more informative output (https://github.com/david-juan/ChromDet/tree/master/S3Det_modified/)
-  Scripts for input pre-processing, S3det execution and output formating (https://github.com/david-juan/ChromDet/tree/master/scripts/)
-  Test files for performing an example toy analysis involving 9 samples from the original article (Carrillo-de Santa Pau et al.) (https://github.com/david-juan/ChromDet/tree/master/test/)

### Steps

#### 1) Download ChromDet (https://github.com/david-juan/ChromDet/archive/master.zip)

#### 2) Uncompress ChromDet
-  ./unzip ChromDet-master.zip

#### 3) Install S3det
-  cd ChromDet/S3Det_modified/
-  make

## Running S3det

### Files

-  Master Script for Sdet analyses (https://github.com/david-juan/ChromDet/tree/master/scripts/run_S3det_analysis.pl)
-  Pre-processing script (https://github.com/david-juan/ChromDet/tree/master/scripts/prepare_S3det_analysis.pl)
-  Script for running S3det and formating its output (https://github.com/david-juan/ChromDet/tree/master/scripts/S3det_interface.pl)
-  Bedfiles containing chromating states segmentations for the samples to analyse (nine examples at https://github.com/david-juan/ChromDet/tree/master/test/).
-  Chromatin states collapses definition (an example at https://github.com/david-juan/ChromDet/tree/master/test/States_collapse.txt)
-  Samples human readable names (an example at https://github.com/david-juan/ChromDet/tree/master/test/Samples_beds.tsv)

#### Output Files

##### Main output files

- Chromatin space file. A tsv file containing the coordinates of the samples in the Chromatin Space, as well as to what cluster tehy belong
- Chromatin Determinants Regions (CDRs) File. A tsv file containing the genomic coordinates of the CDRs, as well as to the split of which clusters they contribute and their pattern of collapsed chromatin states in the different samples.

##### Other output files

- Whole genome patterns of collapsed chromatin states. A tsv file containing the genomic coordinates of every region of the genome characterized by a different pattern of collapsed chromatin states along different samples than their neighbouring regions.
- Filtered regions with patterns of collapsed chromatin states. A tsv containing the genomic coordinates of every region of the genome with a pattern of collapsed chromatin states along different samples considered informative.
- Fasta file of the filtered regions with patterns of collapsed chromatin states. A fasta file, where regions with patterns of collapsed chromatin states along different samples are represented as gap-free aligned sequences of chromatin states for every sample.
- S3det raw resuls. File with the raw output for S3det


### External software requirements

-  Pre-processing script requires bedtools 2 or higher to be installed (http://bedtools.readthedocs.io/en/latest/)

### Steps

#### 1) Go to scripts dir

-  cd ChromDet/scripts

#### 2) Run the complete analysis (minimal)

-  ./run_S3det_analysis.pl -d \<dir containing the bedfiles\>

### 3) Run the complete analysis (recommended)

- ./run_S3det_analysis.pl -d \<dir containing the bedfiles\> -c \<file with the equivalences between chromatin states and their collapse\> -a \<file with human meaningful names for the samples\> -v

#### Toy Example

-  ./run_S3det_analysis.pl -d ../test/ -c ../test/States_collapse.txt -a ../test/Samples_beds.tsv -s ../S3Det_modified/ -v
