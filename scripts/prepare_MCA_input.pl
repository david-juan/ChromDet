#!/opt/local/bin/perl

#This script transforms individual bedfiles containing genome segmentations for different samples in a fasta file and bedfile containing those regions suitable for S3det, a MCA-based approach for analysis of intersample epigenomic variability

use Getopt::Long;
use Pod::Usage;
use Statistics::R;
use strict;
#use warnings;

=head1 NAME

prepare_MCA_input.pl

=head1 SYNOPSIS

    Usage /path/to/prepare_MCA_input.pl [options] 

     Options:
	   -d		'bed_dir' [MANDATORY] directory containing the set of genome segmentations in bed format.
	   -a		'sample annotation file' [OPTIONAL] file containing the list of names of the bed files to analyse
	 			(mnemonic names of the samples to include in the results can be included in a second tab separated column).
	 			In case this file is not provided, all bed file in 'bed_dir' will be used. 
	   -n		'min_num_states' [OPTIONAL] minimum number of states found in a region along all the samples (default value = 2)
	   -m		'min_num_samples_per_state' [OPTIONAL] minimum number of samples with an state to be considered in filtering (see -n; default value = 2)
	   -r		'min_num_regions_pattern' [OPTIONAL] minimal number of regions with the same pattern of states along the samples to be included in the analysis (default value = 10)
	   -o		'output_prefix' [OPTIONAL] prefix added to the output files (default = "collapsed_samples")
	   -b		'bedtools_dir' [OPTIONAL] path to bedtools (by default uses $PATH)
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

#/Users/dadejuan/SMarsili/eggNOG/getting_ortholog_matrix.pl all /Users/dadejuan/SMarsili/eggNOG/eggNOG_metazoa_species_taxid.txt dataB/ ../../ /Users/dadejuan/SMarsili/tree_to_matrix.R /Users/dadejuan/SMarsili/eggNOG/tax_sciname-reformat.txt


my ($bed_dir,$beds_file,$opt_help,$bedtools_path,$command,$cnt,$pre_cnt,$header,$prev_chr,$prev_start,$prev_end,$prev_pattern,$prev_ok,$sample,$state,$collapsed_sample_filtered);
my ($prev_ok_region,$min_states,$min_samples_states,$min_regions_pattern,$out_pre,$bedtools_path);
my (@AA,@tr,@tr2,@bed_files,@working_beds,@sample_names,@prev_ok_regions,@names,@seq,$i,,@states);
my (%working_mnemo_bed_files,%states,%cont_patterns,%all_states,%states2code);

$min_states=2;
$min_samples_states=2;
$min_regions_pattern=10;
$out_pre="samples";
$bedtools_path='';
@AA=('K','V','g','i','j','k','l','m','n','A','B','C','D','E','F','G','H','I','J','L','M','N','O','P','Q','R','S','T','U','W','X','Y','Z','a','b','c','d','e','f','h');
# Get commandline arguments
GetOptions (
			'd=s' => \$bed_dir,
			'a=s' => \$beds_file,
			'n=i' => \$min_states,
			'm=i' => \$min_samples_states,
			'r=i' => \$min_regions_pattern,
			'o=s' => \$out_pre,
			'b=s' => \$bedtools_path,
			'help!' =>  \$opt_help
			) or pod2usage( "Try '$0 --help' for more information." ) && exit;

pod2usage( -verbose => 2 ) if $opt_help || !$bed_dir;

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
opendir BED_DIR, "$bed_dir" or die "I couldn't open $bed_dir directory\n";
while(my $bed_file=readdir(BED_DIR))
{
	if($bed_file =~/\.bed/ && -s "$bed_dir/$bed_file")
	{
		if(!$beds_file or exists($working_mnemo_bed_files{$bed_file}))
		{
			$pre_cnt=$cnt;
			$cnt++;
			if(exists($working_mnemo_bed_files{$bed_file})){push @sample_names, $working_mnemo_bed_files{$bed_file};}
			else{push @sample_names, $bed_file;}
			if(!$command)
			{
				$command="$bedtools_path makewindows -b $bed_dir/$bed_file -w 200 -i src > $bed_dir/tmp_$cnt";
			}else
			{
				$command="$bedtools_path intersect -a $bed_dir/tmp_$pre_cnt -b $bed_dir/$bed_file -wao |cut -f1-4,8 | perl -ne \'s\/\\t(\\w\+\\n\)\/\|\$1\/;print;\'  > $bed_dir/tmp_$cnt;rm $bed_dir/tmp_$pre_cnt";
			}
			`$command`;
		}
	}
}
closedir BED_DIR;


$"="|";
$header="\"#Chr\tStart\tEnd\t@sample_names\"";
push @prev_ok_regions,"#Chr\tStart\tEnd\t@sample_names\n";
system "echo $header > $bed_dir/$out_pre\_collapsed.tab;cat $bed_dir/tmp_$cnt >> $bed_dir/$out_pre\_collapsed.tab";
system "rm $bed_dir/tmp_$cnt";

$prev_ok=0;
open COLLAPSED_FILE, "$bed_dir/$out_pre\_collapsed.tab" or die "I couldn't open $bed_dir/$out_pre\_collapsed.tab\n";
while(<COLLAPSED_FILE>)
{
	chomp;
	@tr=split /\t/;
	if($prev_chr eq $tr[0] && $prev_pattern eq $tr[$#tr])
	{
		if($prev_ok){$prev_end=$tr[1];}
		else{next;}
	}else
	{
		if($prev_ok){push @prev_ok_regions, "$prev_chr\t$prev_start\t$prev_end\t$prev_pattern";}
		$prev_ok=0;
		if($cont_patterns{$tr[3]}>$min_regions_pattern)
		{
			$prev_ok=1;
			$prev_chr=$tr[0];
			$prev_start=$tr[1];
			$prev_end=$tr[2];
			$prev_pattern=$tr[3];
			$cont_patterns{$tr[3]}++;
			next;
		}
		@tr2=split/\|/,$tr[$#tr];
		undef %states;
		foreach $sample(@tr2)
		{
			$states{$sample}++;
		}
		$cnt=0;
		foreach $state(keys(%states))
		{
			if(!exists($all_states{$state})){$all_states{$state}=1;}
			if($states{$state}>$min_samples_states)
			{
				$cnt++;
			}
		}
		if($cnt>$min_states)
		{
			$prev_ok=1;
			$prev_chr=$tr[0];
			$prev_start=$tr[1];
			$prev_end=$tr[2];
			$prev_pattern=$tr[3];
			$cont_patterns{$tr[3]}++;
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
open OUT_BED_FILE, ">$bed_dir/$out_pre\_collapsed_filtered.tab" or die "I couldn't open $bed_dir/$out_pre\_collapsed_filtered.tab\n";;
foreach $prev_ok_region(@prev_ok_regions)
{
	@tr=split/\t/,$prev_ok_region;
	if($prev_ok_region =~ /^#/){print OUT_BED_FILE $prev_ok_region;chomp $tr[$#tr];@names=split/\|/,$tr[$#tr];next;}
	if($cont_patterns{$tr[3]}>$min_regions_pattern)
	{
		print OUT_BED_FILE "$prev_ok_region\t$cont_patterns{$tr[3]}\n";
		@tr=split/\|/,$tr[3];
		for($i=0;$i<@tr;$i++)
		{
			$seq[$i].=$states2code{$tr[$i]};
		}
	}
}
close OUT_BED_FILE;



open OUT_FASTA_FILE, ">$bed_dir/$out_pre\_collapsed_filtered.fa" or die "I couldn't open $bed_dir/$out_pre\_collapsed_filtered.fa\n";
for($i=0;$i<@names;$i++)
{
	print OUT_FASTA_FILE ">$names[$i]\n$seq[$i]\n";
}
close OUT_FASTA_FILE;





