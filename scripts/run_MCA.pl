#!/opt/local/bin/perl

#This script run S3det out of the fasta file obtained by prepare_MCA_input.pl and uses the index file obtained from the same script to parse S3det output file retrieving a file with the samples coordinates in the selected MCA principal components and another file with the Chromatin Determining Regions

use Getopt::Long;
use Pod::Usage;
use Statistics::R;
use strict;
#use warnings;

=head1 NAME

run_MCA.pl

=head1 SYNOPSIS

    Usage /path/to/run_MCA.pl [options] 

     Options:
	   -d		'run_dir' [MANDATORY] directory containing the the fasta and index files.
	   -f		'fasta file' [OPTIONAL] file containing the a one-character codification of the analysed regions obtained by prepare_MCA_input.pl (default: it uses the first fasta file  detected in 'run_dir').
	   -x		'index_file' [OPTIONAL] file containing the index file of the analysed regions obtained by prepare_MCA_input.pl (default: it uses the first index file detected in 'run_dir').
	   -o		'output_prefix' [OPTIONAL] prefix added to the output files (default = "")
	   -s		's3det_path' [OPTIONAL] path to S3det (by default uses $PATH)
	   -c		's3det_opts'  [OPTIONAL] quotes delimited string character with the options used for S3det (default = "-v"). See S3det documentation for information about available options
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


my ($run_dir,$fasta_file,$index_file,$s3det_opts,$out_pre,$s3det_path,$pre_s3det_opts,$opt_help,$file,$num_axes,$i,$general_out_file,$header,$s3det_file,$sdps_file);
#my ();
my (@tr,@coords,@sdps,@index_file,@tr2,@reorder_samples,@new_samples,@tmp);
my (%cluster_samples,%samples,%coords);


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
			'help!' =>  \$opt_help
			) or pod2usage( "Try '$0 --help' for more information." ) && exit;

pod2usage( -verbose => 2 ) if $opt_help || !$run_dir;

if($pre_s3det_opts){$s3det_opts.=" $pre_s3det_opts";}
if($s3det_path)
{
	$s3det_path.="/S3det_v2.2.exe";
}else
{
	$s3det_path="S3det_v2.2.exe";
}


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
	$general_out_file=$fasta_file;
	$general_out_file=~s/\.\w+$//;
	$s3det_file=$general_out_file.".s3";
	if(system "$s3det_path -i $run_dir/$fasta_file -o $run_dir/$s3det_file $s3det_opts")
	{
		open S3DET_FILE, "$run_dir/$s3det_file" or die "I couldn't open $run_dir/$s3det_file\n";
		while(<S3DET_FILE>)
		{
			if(/^UI: Number of axes selected: (\d+)\n/){$num_axes=$1;}
			elsif(/^SeqCoord\:/)
			{
				@tr=split/\s+/;
				$i=$num_axes+1;
				push @coords{$1},@tr[2..$i];
			}elsif(/^CL\:\s+(\S+)\s+(\d+)\n/){push @{$cluster_samples{$2}},$1;}
			elsif(/^RE\:/){@tr=split/\s+/;push @sdps,$tr[1];};
		}
		close S3DET_FILE;
		if(!$num_axes || !@coords || !keys(%cluster_samples)){die "S3det input file is wrong\n"}
		elsif(!@sdps){print STDERR "I didn't detect any CDR\n";}
		open INDEX_FILE, "run_dir/$index_file" or die "I couldn't open run_dir/$fasta_file\n";
		@index_file=<INDEX_FILE>;
		close INDEX_FILE;
		if($index_file[0]=~/^#/)
		{
			chomp $index_file[0];
			@tr=split /\t/,$index_file[0];
			if(@tr!=4){die "ERROR: I couldn't find a correct header in $index_file\n";}
			@tr2=split/\|/;$tr[$#tr];
			if(@tr2 !=  @{$cluster_samples{"1"}}){die "ERROR: I couldn't find a correct header in $index_file\n";}
			for($i=0;$i<@tr2;$i++)
			{
				$samples{$tr2[$i]}=$i;
			}
			for($i=1;@{$cluster_samples{$i}};$i++)
			{
				foreach (@{$cluster_samples{$i}})
				{
					push @reorder_samples,$samples{$_};
					push @new_samples,$_;
				}
			}
			$"="\t";
			$header="@tr[0..@tr]\t";
			$"="|";
			$header.="@new_samples\n";
		}else{die "ERROR: I couldn't find a correct header in $index_file\n";}
		for($i=1;@{$cluster_samples{$i}};$i++)
		{
			foreach (@{$cluster_samples{$i}})
			{
				push @reorder_samples,$samples{$_};
				push @new_samples,$_;
			}
		}
		$sdps_file=$general_out_file.".sdp";
		open SDPS_FILE, ">$run_dir/$sdps_file" or die "I couldn't open $run_dir/$sdps_file\n";
		print SDPS_FILE $header;
		foreach (@sdps)
		{
			@tr=split /\t/, $index_file[$_];
			@tr2=split/\|/;$tr[$#tr];
			$"="\t";
			print SDPS_FILE "@tr[0..@tr]\t";
			undef @tmp;
			foreach (@reorder_samples)
			{
				push @tmp, $tr2[$_];
			}
			$"="|";
			print SDPS_FILE "@tmp\n";
			
		}
		close SDPS_FILE
	}else{die "S3det execution \"$s3det_path -i $run_dir/$fasta_file -o $run_dir/$s3det_file $s3det_opts\" failed\n"};
}else
{
	die "I couldn't find adequate fasta and index files in $run_dir for running S3det\n";
}

