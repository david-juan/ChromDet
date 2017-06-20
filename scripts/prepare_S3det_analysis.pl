#!/usr/bin/perl

#This script transforms individual bedfiles containing genome segmentations for different samples in a fasta file and bedfile containing those regions suitable for S3det, a MCA-based approach for analysis of intersample epigenomic variability

use Getopt::Long;
use Pod::Usage;
use strict;

=head1 NAME

prepare_S3det_analysis.pl

=head1 SYNOPSIS

    Usage /path/to/prepare_S3det_analysis.pl [options] 

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
	   -v		displays messages to STDERR with real time information about pre-processing
	   --help	prints this brief help message

=head1 DESCRIPTION

This script transform individual bedfiles containing genome segmentations for different samples in a fasta file and a bed file containing those regions suitable for S3det, a MCA-based approach for analysis of intersample epigenomic variability
A modified version of S3det should be provided with this script.


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
my ($prev_ok_region,$min_states,$min_samples_states,$min_regions_pattern,$out_pre,$bedtools_path,$collapses_file,$original_state,$traslated_pattern,$pre_traslation,$verbose);
my ($number,$pattern_num,$prev_pattern_num,$autosomal);
my (@AA,@tr,@tr2,@bed_files,@working_beds,@sample_names,@prev_ok_regions,@names,@seq,$i,@states,@numbers);
my (%working_mnemo_bed_files,%states,%cont_patterns,%all_states,%states2code,%states_collapses,%pattern_collapse,%all_collapse_states,%state_number);

$min_states=2;
$min_samples_states=2;
$min_regions_pattern=10;
$out_pre="samples";
$bedtools_path='';
$autosomal='T';
@AA=('K','V','g','i','j','k','l','m','n','A','B','C','D','E','F','G','H','I','J','L','M','N','O','P','Q','R','S','T','U','W','X','Y','Z','a','b','c','d','e','f','h');
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
			'v!' => \$verbose,
			'help!' =>  \$opt_help
			) or pod2usage( "Try '$0 --help' for more information." ) && exit;

pod2usage( -verbose => 2 ) if $opt_help || !$bed_dir;
$autosomal=uc($autosomal);
if($min_states>0)
{
	$min_states--;
}
if($min_samples_states>0)
{
	$min_samples_states--;
}
if($min_regions_pattern>0)
{
	$min_regions_pattern--;
}
if($bedtools_path)
{
	$bedtools_path.="/bedtools";
}else
{
	$bedtools_path="bedtools";
}

if($beds_file)
{
	open BEDS_FILE, "$beds_file" or die "I couldn't open $beds_file\n";
	while(<BEDS_FILE>)
	{
		chomp;
		@tr=split/\t/;
		push @working_beds, $tr[0];
		if($tr[1]){$working_mnemo_bed_files{$tr[0]}=$tr[1];}
		else{$working_mnemo_bed_files{$tr[0]}=$tr[0];}
	}
	close BEDS_FILE;
}
$cnt=0;
if($verbose){print STDERR "PREPARING - Step 1: Merging bedfiles. Be patient, this step can last a while depending on the number of bedfiles\n";}
opendir BED_DIR, "$bed_dir" or die "I couldn't open $bed_dir directory\n";
while(my $bed_file=readdir(BED_DIR))
{
	if($bed_file =~/\.bed/ && -s "$bed_dir/$bed_file")
	{
		if(!$beds_file or exists($working_mnemo_bed_files{$bed_file}))
		{
			$pre_cnt=$cnt;
			$cnt++;
			if($verbose){print STDERR "\tProcessing file $cnt: $bed_file\n";}
			if(exists($working_mnemo_bed_files{$bed_file})){push @sample_names, $working_mnemo_bed_files{$bed_file};}
			else{push @sample_names, $bed_file;}
			if(!$command)
			{
				$command="$bedtools_path makewindows -b $bed_dir/$bed_file -w 200 -i src > $bed_dir/tmp_$cnt";
			}else
			{
				$command="$bedtools_path intersect -a $bed_dir/tmp_$pre_cnt -b $bed_dir/$bed_file -wao |cut -f1-4,8 | perl -ne \'s\/\\t(\\w\+\\n\)\/\|\$1\/;print;\'  > $bed_dir/tmp_$cnt;rm $bed_dir/tmp_$pre_cnt";
			}
			system "$command";
		}
	}
}
closedir BED_DIR;

if($collapses_file)
{
	open COLLAPSES_FILE, "$collapses_file" or die "I couldn't open $collapses_file\n";
	while(<COLLAPSES_FILE>)
	{
		chomp;
		@tr=split/\t/;
		if(exists($states_collapses{$tr[0]}) and $states_collapses{$tr[0]} ne $tr[1])
		{
			die "ERROR: $collapses_file file contains abiguous collapses: $tr[0] is assigned to $states_collapses{$tr[0]} and to $tr[1]\n";
		}
		$states_collapses{$tr[0]}=$tr[1];
		$all_collapse_states{$tr[1]}=1;
	}
	close COLLAPSES_FILE;
	foreach (keys(%all_collapse_states))
	{
		if(exists($states_collapses{$_}))
		{
			die "ERROR: $collapses_file file contains a conflicting definition of state $_ as input and as output chromatin states.\n";
		}
	}
}


$"="|";
$header="\"#Chr\tStart\tEnd\t@sample_names\t@sample_names\"";
push @prev_ok_regions,"#Chr\tStart\tEnd\t@sample_names\t@sample_names\n";
system "echo $header > $bed_dir/$out_pre\_collapsed.tab;cat $bed_dir/tmp_$cnt >> $bed_dir/$out_pre\_collapsed.tab";
system "rm $bed_dir/tmp_$cnt";

if($verbose){print STDERR "PREPARING - Step 2: Filtering uninformative regions\n";}

$cnt=0;
$prev_ok=0;
open COLLAPSED_FILE, "$bed_dir/$out_pre\_collapsed.tab" or die "I couldn't open $bed_dir/$out_pre\_collapsed.tab\n";
while(<COLLAPSED_FILE>)
{
	chomp;
	@tr=split /\t/;
	if($autosomal == 'T' && !/^\#/ && $tr[0] !~ /chr\d\d*$/){next;}
	if($collapses_file)
	{
		if(!/^\#/ && !exists($pattern_collapse{$tr[3]}))
		{
			$pre_traslation=$tr[3];
			foreach $original_state(keys(%states_collapses))
			{
				$tr[3]=~s/^$original_state\|/$states_collapses{$original_state}|/ig;
				$tr[3]=~s/\|$original_state$/|$states_collapses{$original_state}/ig;
				$tr[3]=~s/\|$original_state\|/|$states_collapses{$original_state}|/ig;
				$tr[3]=~s/\|$original_state\|/|$states_collapses{$original_state}|/ig;
				
			}
			$pattern_collapse{$pre_traslation}=$tr[3];

		}else
		{
			$tr[3]=$pattern_collapse{$tr[3]};
		}
	}
	if($prev_chr eq $tr[0] && $prev_pattern eq $tr[$#tr])
	{
		if($prev_ok){$prev_end=$tr[1];}
		else{next;}
	}else
	{
		@tr2=split/\|/,$tr[3];
		undef %state_number;
		$number=0;
		undef @numbers;
		undef %states;
		foreach (@tr2)
		{
			if(!exists($state_number{$_})){$number++;$state_number{$_}=$number;}
			$states{$_}++;
			push @numbers, $state_number{$_};
		}
		$"="|";
		$pattern_num="@numbers";
		if($prev_ok)
		{
			push @prev_ok_regions, "$prev_chr\t$prev_start\t$prev_end\t$prev_pattern\t$prev_pattern_num";
		}
		$prev_ok=0;
		if($cont_patterns{$pattern_num}>$min_regions_pattern)
		{
			$prev_ok=1;
			$prev_chr=$tr[0];
			$prev_start=$tr[1];
			$prev_end=$tr[2];
			$prev_pattern=$tr[3];
			$prev_pattern_num=$pattern_num;
			$cont_patterns{$pattern_num}++;
			next;
		}
		if(!/\#/)
		{
			$cnt=0;
			foreach $state(keys(%states))
			{
				if(!exists($all_states{$state}))
				{
					$all_states{$state}=1;
				}
				if($states{$state}>$min_samples_states)
				{
					$cnt++;
				}
			}
		}
		if($cnt>$min_states)
		{
			$prev_ok=1;
			$prev_chr=$tr[0];
			$prev_start=$tr[1];
			$prev_end=$tr[2];
			$prev_pattern=$tr[3];
			$prev_pattern_num=$pattern_num;
			$cont_patterns{$pattern_num}++;
		}else
		{
			$prev_ok=0;
		}		
	}
}
close COLLAPSED_FILE;

@states=sort keys(%all_states);
for($i=0;$i<@states;$i++)
{
	$states2code{$states[$i]}=$AA[$i];
}

if($verbose){print STDERR "PREPARING - Step 3: Writing genomic coordinates file ($bed_dir/$out_pre\_collapsed_filtered.tab) for informative regions\n";}

open OUT_BED_FILE, ">$bed_dir/$out_pre\_collapsed_filtered.tab" or die "I couldn't open $bed_dir/$out_pre\_collapsed_filtered.tab\n";;
foreach $prev_ok_region(@prev_ok_regions)
{
	@tr=split/\t/,$prev_ok_region;
	if($prev_ok_region =~ /^#/)
	{
		print OUT_BED_FILE $prev_ok_region;
		chomp $tr[$#tr];
		@names=split/\|/,$tr[$#tr];
		next;
	}
	if($cont_patterns{$tr[4]}>$min_regions_pattern)
	{
		print OUT_BED_FILE "$prev_ok_region\n";
		@tr=split/\|/,$tr[3];
		for($i=0;$i<@tr;$i++)
		{
			$seq[$i].=$states2code{$tr[$i]};
		}
	}
}
close OUT_BED_FILE;


if($verbose){print STDERR "PREPARING - Step 4: Writing fasta file ($bed_dir/$out_pre\_collapsed_filtered.fa) for feeding S3det\n";}
open OUT_FASTA_FILE, ">$bed_dir/$out_pre\_collapsed_filtered.fa" or die "I couldn't open $bed_dir/$out_pre\_collapsed_filtered.fa\n";
for($i=0;$i<@names;$i++)
{
	print OUT_FASTA_FILE ">$names[$i]\n$seq[$i]\n";
}
close OUT_FASTA_FILE;





