#!/usr/bin/perl

#This script run an standard epigenomic S3det analysis. S3det is a MCA-based approach for analysis of intersample epigenomic variability. The expected input is a set of bed files containing chromatin states segmentations for different samples from different cell types of the same species. A number of options are available to help designing an experiment-specific analysis and to provide easy to interpret results. The expected output is a Chromatin Sample Space, reflecting the sample clustering and a list of Chromatin Determining Positions (CDRs) that are the major contributors to this sample clustering and therefore are expected to be important chromatin regions for defining the epigenomic framework characterizing every cluster in oposition to all the rest of the clusters.

use Getopt::Long;
use Pod::Usage;
use strict;
use Cwd 'abs_path';

=head1 NAME

run_S3det_analysis.pl

=head1 SYNOPSIS

    Usage /path/to/run_S3det_analysis.pl [options] 

     Options:
	   -d		'bed_dir' [MANDATORY] directory containing the set of genome segmentations in bed format.
	   -a		'sample annotation file' [OPTIONAL] file containing the list of names of the bed files to analyse
	 			(mnemonic names of the samples to include in the results can be included in a second tab separated column).
	 			In case this file is not provided, all bed file in 'bed_dir' will be used.
	   -c		'state_collapses_file' [OPTIONAL] two-column tabular file containg equivalences between each chromatin state present in the genome segmentations and collapsed states.
	 			Collapses are recommended to be based on their shared biological role (eg. different states replecting enhancers).
	   -h		'T/F' [OPTIONAL] 'T' for selecting only autosomal chromosomes. 'F' for using all chromosomes (default=T)
	   -n		'min_num_states' [OPTIONAL] minimum number of states found in a region along all the samples (default value = 2)
	   -m		'min_num_samples_per_state' [OPTIONAL] minimum number of samples with an state to be considered in filtering (see -n; default value = 2)
	   -r		'min_num_regions_pattern' [OPTIONAL] minimal number of regions with the same pattern of states along the samples to be included in the analysis (default value = 10)
	   -o		'output_prefix' [OPTIONAL] prefix added to the output files (default = "collapsed_samples")
	   -b		'bedtools_dir' [OPTIONAL] path to bedtools (by default uses $PATH)
	   -s		's3det_path' [OPTIONAL] path to S3det (by default uses $PATH)
	   -f		's3det_opts'  [OPTIONAL] quotes delimited string character with the options used for S3det (default = "-v"). See S3det documentation for information about available options
	   -v		displays messages to STDERR with real time information about the ongoing analysis
	   --help	prints this brief help message

=head1 DESCRIPTION

This script run an standard epigenomic S3det analysis. S3det is a MCA-based approach for analysis of intersample epigenomic variability. The expected input is a set of bed files containing chromatin states segmentations for different samples from different cell types of the same species. A number of options are available to help designing an experiment-specific analysis and to provide easy to interpret results. The expected output is a Chromatin Sample Space, reflecting the sample clustering and a list of Chromatin Determining Positions (CDRs) that are the major contributors to this sample clustering and therefore are expected to be important chromatin regions for defining the epigenomic framework characterizing every cluster in oposition to all the rest of the clusters.

=head1 LICENSE

This is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
    
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details (<http://www.gnu.org/licenses/>).

=head1 AUTHOR

David Juan (david.juan@upf.edu)>

=cut




$"="\t";

my ($bed_dir,$beds_file,$opt_help,$bedtools_path,$command,$cnt,$pre_cnt,$header,$prev_chr,$prev_start,$prev_end,$prev_pattern,$prev_ok,$sample,$state,$collapsed_sample_filtered);
my ($prev_ok_region,$min_states,$min_samples_states,$min_regions_pattern,$out_pre,$bedtools_path,$collapses_file,$original_state,$traslated_pattern,$current_dir);
my ($pre_s3det_opts,$s3det_path,$verbose,$pre_command,$abs_bed_dir,$run_command,$autosomal);
my (@AA,@tr,@tr2,@bed_files,@working_beds,@sample_names,@prev_ok_regions,@names,@seq,$i,@states);
my (%working_mnemo_bed_files,%states,%cont_patterns,%all_states,%states2code,%states_collapses,%pattern_collapse,%all_collapse_states);

#$min_states=2;
#$min_samples_states=2;
#$min_regions_pattern=10;
$out_pre="samples";
#$bedtools_path='';
$s3det_path="./";
$autosomal='T';
# Get commandline arguments
GetOptions (
			'd=s' => \$bed_dir,
			'a=s' => \$beds_file,
			'c=s' => \$collapses_file,
			'h=s' => \$autosomal,
			'n=i' => \$min_states,
			'm=i' => \$min_samples_states,
			'r=i' => \$min_regions_pattern,
			'o=s' => \$out_pre,
			'b=s' => \$bedtools_path,
			'f=s' => \$pre_s3det_opts,
			's=s' => \$s3det_path,
			'v!' => \$verbose,
			'help!' =>  \$opt_help
			) or pod2usage( "Try '$0 --help' for more information." ) && exit;

pod2usage( -verbose => 2 ) if $opt_help || !$bed_dir;

$autosomal=uc($autosomal);
$pre_command = "./prepare_S3det_analysis.pl -d $bed_dir";
if($beds_file){$pre_command .=" -a $beds_file";}
if($collapses_file){$pre_command .=" -c $collapses_file";}
if($min_states){$pre_command .=" -n $min_states";}
if($min_samples_states){$pre_command .=" -m $min_samples_states";}
if($min_regions_pattern){$pre_command .=" -r $min_regions_pattern";}
if($out_pre){$pre_command .=" -o $out_pre";}
if($bedtools_path){$pre_command .=" -b $bedtools_path";}
if($autosomal != 'T'){$pre_command .=" -h $autosomal";}
if($verbose){$pre_command .=" -v";}

if($verbose){print STDERR "Running $pre_command\n";}
system qq[$pre_command];

$current_dir = abs_path();
$current_dir =~ s/\s/\\\ /ig;
$abs_bed_dir = abs_path($bed_dir);
$abs_bed_dir =~ s/\s/\\\ /ig;
if($pre_s3det_opts)
{
	$run_command = qq [$current_dir/S3det_interface.pl -d $abs_bed_dir -f $out_pre\_collapsed_filtered.fa -x $out_pre\_collapsed_filtered.tab -s $s3det_path -c "$pre_s3det_opts"];
}else
{	
	$run_command = qq [$current_dir/S3det_interface.pl -d $abs_bed_dir -f $out_pre\_collapsed_filtered.fa -x $out_pre\_collapsed_filtered.tab -s $s3det_path];
}

if($verbose)
{
	$run_command.= " -v";
	print STDERR "\nRunning $run_command\n";
}

system "$run_command";


