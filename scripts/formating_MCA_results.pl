#!/opt/local/bin/perl

#This script will extract Chromatin Space coordinates for the analysed samples to '<output_prefix>_chromatin_space.tsv' and Chromatin Determining Regions (CDRs) to '<output_prefix_CDRs.tsv' from the S3det results file


use Getopt::Long;
use Pod::Usage;
use Statistics::R;
use strict;
#use warnings;

=head1 NAME

formating_MCA_results.pl

=head1 SYNOPSIS

    Usage /path/to/formating_MCA_results.pl [options] 

     Options:
	   -r		'S3det_results_file' [MANDATORY] file containing the results from S3Det analysis
	   -b		'target_regions_bedfile' [MANDATORY] bedfile containing the regions analysed by S3det. IMPORTANT: In order to ensure coherence to the input fasta file for S3Det, both files should be generated together by prepare_MCA_input.pl and they must remain unmodified.
	   -o		'output_prefix' [OPTIONAL] prefix added to the output files (default = "S3det_results_file without file extension")
	   --help	prints this brief help message

=head1 DESCRIPTION

This script will extract Chromatin Space coordinates for the analysed samples to '<output_prefix>_chromatin_space.tsv' and Chromatin Determining Regions (CDRs) to '<output_prefix_CDRs.tsv' from the S3det results file
A modified version of S3det suitable for epigenomic analyses should be provided with this script.


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



my ($i,$results_file,$bed_file,$opt_help,$axes_number,$clusters_number,$out_space,$out_CDRs,$line);
my (@tr,@tr2,@bed_file_lines,@idx_out_samples,@chrom_space,@out_samples,@CDRs);
my (%variance,%seq_cluster,%CDR_split);

# Get commandline arguments
GetOptions (
			'r=i' => \$results_file,
			'o=s' => \$out_pre,
			'b=s' => \$bed_file,
			'help!' =>  \$opt_help
			) or pod2usage( "Try '$0 --help' for more information." ) && exit;

pod2usage( -verbose => 2 ) if $opt_help || !$results_file || $bed_file;

if(!$out_pre)
{
	$out_pre=$results_file;
	$out_pre=~s/\.[\w]$//;
}

$out_space=$out_pre."_chromatin_space.tsv";
$out_CDRs=$out_pre."_CDRs.tsv";

open BED_FILE, "$bed_file" or die "I couldn't open $bed_file\n";
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

open RESULTS_FILE, "$results_file" or die "I couldn't open $results_file\n";
while(<RESULTS_FILE>)
{
	if(/SeqCoord:/)
	{
		chomp;
		@tr=split/\s+/;
		$"="\t";
		push @chrom_space, "$tr[1..$#tr]\n"; 
	}
	if(/UI\: Number of axes selected: (\d+)\n/){$axes_number=$1;}
	if(/UI: Percentage of variance explained by the axis (\d+)  \(selected\): ([\d\.]+) %/){$variance{$1}=$2;}
	if(/SC: NO/){die "No stable clustering was detected\n";}
	if(/UI: Number of groups selected: (\d+)\n/){$clusters_number=$1;}
	elsif(/^CL\s+(\S+)\s+(\d+)\n/){$seq_cluster{$1}=$2;}
	elsif(/^RE:/)
	{
		@tr=split/\s+/;
		$"="|";
		$CDR_split{$tr[1]}="@tr[5..$#tr]";
	}
	
}
close RESULTS_FILE;

open SPACE_OUT, ">$out_space" or die "I couldn't open $out_space\n";

print SPACE_OUT "Sample";
for($i=1;$i<=$axes_number;$i++)
{
	print SPACE_OUT "\tPC$i ($variance{$i}%)";
}
print SPACE_OUT "\tCluster\n";

foreach(@chrom_space)
{
	@tr=split /\t/;
	print SPACE_OUT "$tr[0]";
	for($i=1;$i<=$axes_number;$i++)
	{
		print SPACE_OUT "\t$tr[$i]";
	}
	print SPACE_OUT "Cluster_$seq_cluster{$tr[0]}\n";
}
close SPACE_OUT;

open CDRS_OUT, ">$out_CDRS" or die "I couldn't open $out_CDRs\n";

@idx_out_samples=sort{$seq_cluster{$bed_samples[$a]}<=>$seq_cluster{$bed_samples[$b]}}0 .. $#bed_samples;
@out_samples=@bed_samples[@idx_out_samples];

print CDRS_OUT "#Chr\tStart\tEnd\tClustersSplit\n";
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

