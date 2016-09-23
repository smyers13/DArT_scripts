#!/usr/bin/perl -w

# dart_filter.pl
# Terry Bertozzi 12.ii.2015
#
# dart_filter_v2.pl
# Steven Myers 4.vii.2016 -- modified to be compatable with recent changes to DArT Report format
#
# Filter DArT RAD data, find duplicates and check the Hamming Distances to identify
# potential problems with the data.
#
# usage: dart_filter.pl --help 
#
# This script is referenced in the publication:
#
# Donnellan, S.C., Foster, R., Junge, C., Huveneers, C., Kilian, A. and Bertozzi, T. (2015) 
# Fiddling with the proof: the Magpie Fiddler Ray is a colour pattern variant of the common 
# Southern Fiddler Ray (Rhinobatidae: Trygonorrhina). _Zootaxa_ (submitted)

# input: .csv output of R script 'DArT_cloneFilter.R'; then potentially .csv output of R script 'DArT_hamFilter.R'

use strict;
use File::Basename;
use Getopt::Long;

my %HoH;
my %sequences; #sequences to test
my %sequences_copy;
my $total_seq=0; # total number of sequences processsed
my ($filt_call, $filt_clust, $filt_adist)=0; #keep track of how many loci are filtered
my $header;
my @duplicates; #clones represented more than once

# getopt variables
my $help;
my $infile; # Name of input file
my $clust; # ClusterSize - numeric
my $adist; # AlleleSeqDist - numeric
my $call; #CallRate - proportion

# declare the command line flags/options
my $opt=GetOptions("adist=i"=>\$adist,
                   "clust=i"=>\$clust,
                   "call=f"=>\$call,
                   "in=s"=>\$infile,
                   'help'=>\$help); 


if (!($opt && $infile && $call) || $help) {#check for the required inputs
    print_help();
}

# check $call is a proportion
if ($call) {
    unless ($call>=0 && $call <=1){
        print_help();
    }
}

unless (open (IN, '<:crlf',$infile)) { #use PerlIO layer to handle windows files
    die ("Can't open input file $infile\n");
}

my ($fname, $dir, $ext) = fileparse($infile,'\..*');
my $outfile_filt=$fname."_filtered.csv";
my $outfile_ham=$fname.".ham_dist_values";
my $outfile_stats=$fname.".stats";

unless(open(OUT,">",$outfile_filt)){
        die ("Unable to open $outfile_filt\n")
}

unless(open(OUT1,">",$outfile_stats)){
        die ("Unable to open $outfile_stats\n")
}

# filter the data based on input options
while (my $line = <IN>) {
    next if ($line=~ /^\*/); #Ignore some header lines
    if ($line=~/^C/){
        $header = $line;
        print OUT $header;
        next;
    }
   
   my ($cloneID, $seq, $snp, $snp_pos, $call_rate, $rest) = split(/,/,$line,6);
    $total_seq++;
    
#    if ($clust) {
#        unless ($clust_size <= $clust){
#            $filt_clust++;
#            next;
#        }
#    }
#    
#    if ($adist) {
#        unless ($allele_dist <= $adist){
#            $filt_adist++;
#            next;
#        }
#    }
#    
    if ($call) {
        unless ($call_rate >= $call){
            $filt_call++;
            next;
        }
    }
    
    print OUT $line;
}

# calculate and print stats
my $running_total = $total_seq/2;
print OUT1 "Total no of loci: ".$running_total;

if($clust){
    print OUT1 "\nLoci removed due to cluster size ($clust): ".$filt_clust/2;
    $running_total=$running_total - $filt_clust/2;
}

if($adist){
    print OUT1 "\nLoci removed due to allele distribution ($adist): ".$filt_adist/2;
     $running_total=$running_total - $filt_adist/2;
}

if($call){
    print OUT1 "\nLoci removed due to call rate ($call): ".$filt_call/2;
    $running_total=$running_total - $filt_call/2;
}

print OUT1 "\nTotal number of loci retained: $running_total\n";
    
close IN;
close OUT;
close OUT1;

unless (open (IN, '<:crlf',$outfile_filt)) { #use PerlIO layer to handle windows files
    die ("Can't open input file $outfile_filt\n");
}

unless(open(OUT,">",$outfile_ham)){
        die ("Unable to open $outfile_ham\n")
}

# check the hamming distances between retained loci
while (my $line = <IN>) {
#edit the following line to match data file (don't forget to change number)
   my ($cloneID, $seq, $snp, $snp_pos, $call_rate, $rest) = split(/,/,$line,6);
    next unless ($snp); #Ignore the reference allele lines with no snp position information
 
    # get a representative of each clone   
    if (exists($sequences{$cloneID})) {
        push @duplicates, $cloneID;
    }else{
        $sequences{$cloneID}=$seq;
    }
}

%sequences_copy=%sequences; #make a copy of the hash

# compare each sequence against all others
foreach my $key (keys %sequences){
    my @result;
    foreach my $key2 (keys %sequences_copy){
        next if($key eq $key2);
        my $result = hamdist($sequences{$key},$sequences_copy{$key2});
        if ($result < 10) {
            $HoH{$key}{$key2} = $result;
        }
    }

    delete $sequences{$key};
    delete $sequences_copy{$key};
}

foreach my $key (sort keys %HoH){
    print OUT "$key: ";
    foreach my $match (sort keys %{$HoH{$key}}){
        print OUT "$match (".$HoH{$key}{$match}.") "
    }
    print OUT "\n";
}

if (@duplicates) {
    print OUT "The following ".scalar(@duplicates)." cloneIDs seem to be duplicated: \n";
    print OUT join(',', sort {$a <=> $b} @duplicates);
}

close IN;
close OUT;

 
# subroutines
sub hamdist{
    return ($_[0] ^ $_[1]) =~ tr/\001-\255//;
}

sub print_help{
    print STDERR "\nExample usage:\n\n";
    print STDERR "$0 --[filters] --in <filename> \n\n";
    print STDERR "      --help = This usage example\n";
    print STDERR "      --in = name of the DArT datafile to filter\n\n";
    print STDERR "   filters:\n";
    print STDERR "      --adist <integer> = sequence distance in cluster (AlleleSeqDist)\n";
    print STDERR "      --clust <integer> = cluster size (CLusterSize)\n";
    print STDERR "      --call <proportion> = minimum coverage across samples (CallRate)\n\n";
    exit;
}
