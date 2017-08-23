#!/usr/bin/perl

#This script run S3det out of the fasta file obtained by prepare_MCA_input.pl and uses the index file obtained from the same script to parse S3det output file retrieving a file with the samples coordinates in the selected MCA principal components and another file with the Chromatin Determining Regions

use Getopt::Long;
use Pod::Usage;

=head1 NAME

S3det_interface.pl

=head1 SYNOPSIS

    Usage /path/to/run_MCA.pl [options] 

     Options:
	   -d		'run_dir' [MANDATORY] directory containing the fasta and index files.
	   -f		'fasta file' [OPTIONAL] file containing the a one-character codification of the analysed regions obtained by prepare_MCA_input.pl (default: it uses the first fasta file  detected in 'run_dir').
	   -x		'index_file' [OPTIONAL] file containing the index file of the analysed regions obtained by prepare_MCA_input.pl (default: it uses the first index file detected in 'run_dir').
	   -o		'output_prefix' [OPTIONAL] prefix added to the output files (default = "")
	   -s		's3det_path' [OPTIONAL] path to S3det (by default uses $PATH)
	   -c		's3det_opts'  [OPTIONAL] quotes delimited string character with the options used for S3det (default = "-v"). See S3det documentation for information about available options
	   -v		displays messages to STDERR with real time information about the ongoing analysis
	   --help	prints this brief help message

=head1 DESCRIPTION

This script run S3det out of the fasta file obtained by prepare_MCA.pl and uses the index file obtained from the same script to parse S3det output file retrieving a file with the samples coordinates in the selected MCA principal components and another file with the Chromatin Determining Regions
prepare_MCA_input.pl script and modified version of S3det should be provided with this script.


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


my ($run_dir,$fasta_file,$index_file,$s3det_opts,$out_pre,$s3det_path,$pre_s3det_opts,$opt_help,$file,$axes_number,$clusters_number,$i,$general_out_file,$s3det_file,);
my ($results_file,$out_bed_line,$out_space,$out_CDRs,$line);
my (@tr,@tr2,@reorder_samples,@new_samples,@tmp,@bed_samples,@bed_file_lines,@idx_out_samples,@chrom_space,@out_samples,@CDRs);
my (%variance,%seq_cluster,%CDR_split,%space_out);


$s3det_path='';
$s3det_opts='-v';
# Get commandline arguments
GetOptions (
			'd=s' => \$run_dir,
			'f=s' => \$fasta_file,
			'x=s' => \$index_file,
			'c=s' => \$pre_s3det_opts,
			'o=s' => \$out_pre,
			's=s' => \$s3det_path,
			'v!' => \$verbose,
			'help!' =>  \$opt_help
			) or pod2usage( "Try '$0 --help' for more information." ) && exit;

pod2usage( -verbose => 2 ) if $opt_help || !$run_dir;

$s3det_path=~ s/\s/\\\ /ig;
$run_dir=~ s/\s/\\\ /ig;
$fasta_file=~ s/\s/\\\ /ig;
$index_file=~ s/\s/\\\ /ig;
$out_pre=~ s/\s/\\\ /ig;
if($pre_s3det_opts){$s3det_opts.=" $pre_s3det_opts";}

if(!$fasta_file || !$index_file)
{
	opendir RUN_DIR, "$run_dir" or die "I couldn't open $run_dir directory\n";
	while(my $file=readdir(RUN_DIR))
	{
		if($file =~/\_collapsed_filtered\.fa/ && -s "$run_dir/$file")
		{
			$fasta_file=$file;
		}elsif($file =~/\_collapsed_filtered\.tab/ && -s "$run_dir/$file")
		{
			$index_file=$file;
		}
	}
	closedir RUN_DIR;
}

if($fasta_file && $index_file)
{
	if(!$out_pre)
	{
		$general_out_file=$fasta_file;
		$general_out_file=~s/\.\w+$//;
		$s3det_file=$general_out_file.".s3";
	}else
	{
		$s3det_file=$out_pre.".s3";
	}
	if($verbose){print STDERR "RUNNING - Step 1: Running $s3det_path/S3det_v2.4.exe -i $run_dir/$fasta_file -o $run_dir/$s3det_file $s3det_opts\n";}
	if(!system "cd $s3det_path;./S3det_v2.4.exe -i $run_dir/$fasta_file -o $run_dir/$s3det_file $s3det_opts")
	{
	
		if($verbose){print STDERR "RUNNING - Step 2: Extracting S3det results\n";}
	
		if(!$out_pre)
		{
			$out_pre=$s3det_file;
			$out_pre=~s/\.\w+$//;
		}

		$out_space=$out_pre."_chromatin_space.tsv";
		$out_CDRs=$out_pre."_CDRs.tsv";

		open BED_FILE, "$run_dir/$index_file" or die "I couldn't open $run_dir/$index_file\n";
		while(<BED_FILE>)
		{
			chomp;
			@tr=split/\t/;
			if(!@bed_file_lines)
			{
				@bed_samples=split/\|/,$tr[$#tr];
			}
			push @bed_file_lines, $_;
		}
		close BED_FILE;

		open RESULTS_FILE, "$run_dir/$s3det_file" or die "I couldn't open $run_dir/$s3det_file\n";
		while(<RESULTS_FILE>)
		{
			if(/SeqCoord:/)
			{
				chomp;
				@tr=split/\s+/;
				$"="\t";
				push @chrom_space, "@tr[1..$#tr]\n"; 
			}
			if(/UI\: Number of axes selected: (\d+)\n/){$axes_number=$1;}
			if(/UI: Percentage of variance explained by the axis (\d+)  \(selected\): ([\d\.]+) %/){$variance{$1}=$2;}
			if(/SC: NO/){die "No stable clustering was detected\n";}
			if(/UI: Number of groups selected: (\d+)\n/){$clusters_number=$1;}
			elsif(/^CL\:\s+(\S+)\s+(\d+)\n/){$seq_cluster{$1}=$2;}
			elsif(/^RE:/)
			{
				s/\(rank\:\d+\)//ig;
				@tr=split/\s+/;
				$"="|";
				$CDR_split{$tr[1]}="@tr[5..$#tr]";
				$CDR_split{$tr[1]}=~s/\_$//;
			}
	
		}
		close RESULTS_FILE;

		open SPACE_OUT, ">$run_dir/$out_space" or die "I couldn't open $run_dir/$out_space\n";

		print SPACE_OUT "#Sample";
		for($i=1;$i<=$axes_number;$i++)
		{
			print SPACE_OUT "\tPC$i ($variance{$i}%)";
		}
		print SPACE_OUT "\tCluster\n";

		foreach(@chrom_space)
		{
			@tr=split /\t/;
			$space_out{$tr[0]}="$tr[0]";
			for($i=1;$i<=$axes_number;$i++)
			{
				$space_out{$tr[0]}.="\t$tr[$i]";
			}
			$space_out{$tr[0]}.= "\tCluster_$seq_cluster{$tr[0]}\n";
		}
		@idx_out_samples=sort{$seq_cluster{$bed_samples[$a]}<=>$seq_cluster{$bed_samples[$b]}}0 .. $#bed_samples;
		@out_samples=@bed_samples[@idx_out_samples];
		
		foreach (@out_samples)
		{
			print SPACE_OUT $space_out{$_};
		}
		close SPACE_OUT;
		
		
		open CDRS_OUT, ">$run_dir/$out_CDRs" or die "I couldn't open $run_dir/$out_CDRs\n";

		print CDRS_OUT "#Chr\tStart\tEnd\tClustersSplit\t@out_samples\n";
		@CDRs=sort{$a<=>$b}keys(%CDR_split);
		foreach $line(@CDRs)
		{
			$out_bed_line=$bed_file_lines[$line];
			@tr=split/\t/,$out_bed_line;
			$"="\t";
			print CDRS_OUT "@tr[0..$#tr-1]";
			@tr2=split/\|/,$tr[$#tr];
			@tr2=@tr2[@idx_out_samples];
			$"="|";
			print CDRS_OUT "\t$CDR_split{$line}\t@tr2\n";
		}
		close CDRS_OUT;
		

	}else{die "S3det execution \"$s3det_path -i $run_dir/$fasta_file -o $run_dir/$s3det_file $s3det_opts\" failed\n"};
}else
{
	die "I couldn't find adequate fasta and index files in $run_dir for running S3det\n";
}

