#!/usr/bin/perl

use strict;
use warnings;

use Bio::SearchIO;
use File::Spec;
use Data::Dumper;

#make sure only one input file is specified
my $len = @ARGV;
if ($len != 3) {
    die "Please specify an input file, a run name, and cutoff e-value!\n";
}
my $infile = shift @ARGV;
my $title = shift @ARGV;
my $cutoff = shift @ARGV;

if (!-f $infile) {
    die "Supplied filepath is not valid: $infile";
}

sub help {
  print "

usage:

match_mites_ggsearch_both-ori.pl  <ggsearch_output_file> <run_name> <E-value_cutoff>

example:
perl match_mites_ggsearch_both-ori.pl aedes_bad-tirs-to-unique.ggsearch_out aedes_bad-to-unique 1e-30

";
  exit 1;
}


my ($volume, $in_dir, $filename) = File::Spec->splitpath($infile);
my $match_out_path = File::Spec->catpath($volume, $in_dir, $title . "_matches.txt");
my $unique_out_path = File::Spec->catpath($volume, $in_dir, $title . "_unique.txt");
open(my $match_out, ">", $match_out_path);
open(my $unique_out, ">", $unique_out_path);

print "Opening  search result file . . . ";
my $fa_aln_obj = Bio::SearchIO->new(-format => 'fasta', -file   => $infile);
print "Opened!\n";
print Dumper($fa_aln_obj);
print "\n\n\n";
while( my $result = $fa_aln_obj->next_result ) {
    #print Dumper($result);
    my $query_id = $result->query_name();
    print "Query:  $query_id\n";
    if( my $hit = $result->next_hit ) {
        my $hsp = $hit->next_hsp;
        my $sig;
        print "Has a hit . . .\t";
        my $hit_id = $hit->name();
        $sig = $hit->significance();
        if (($hsp->percent_identity >= 67.0) and $sig < $cutoff) {
            print "Hit is good match\n";
            print $match_out "$query_id\t$hit_id\t$sig\n";
        }
        else {
            print "Hit is not a good match\n";
            print $unique_out "$query_id\n";
        }
    }
    else {
        print "No hits\n";
        print $unique_out "$query_id\n";
    }
}

close($unique_out);
close($match_out);
