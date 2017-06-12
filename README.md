# ChromDet
A method to classify different cell types from their epigenomic profiles and detect genomic regions determining this classification
This repository contains a slightly modified version of S3det (Rausell et al. PNAS 2010) adapted to perform epigenomic analyses of sample segregation into different cell types. It contains several scripts aimed to facilitate this analysis by providing an interface with S3det and retrieving understable output files containing the Chromatin Sample Space and the Chromatin Determining Regions (CDRs).

## Instalation

### Directories

-  S3det distribution slightly modified to allow a higher number of variables and provide a more informative output (https://github.com/david-juan/ChromDet/tree/master/S3Det_modified/)
-  Scripts for input pre-processing, S3det execution and output formating (https://github.com/david-juan/ChromDet/tree/master/scripts/)
-  Test files for performing an example toy analysis involving 9 samples from the original article (Carrillo-de Santa Pau et al.) (https://github.com/david-juan/ChromDet/tree/master/test/)
