#!/usr/bin/perl

use strict;
use warnings;

use Bio::SearchIO;
use File::Spec;
use Data::Dumper;

#make sure only one input file is specified
my $len = @ARGV;
if ($len != 4) {
    die "Please specify an input file, a run name, cutoff e-value, and search type!\n";
}
my $infile = shift @ARGV;
my $title = shift @ARGV;
my $cutoff = shift @ARGV;
my $type = shift @ARGV;

if (!-f $infile) {
    die "Supplied filepath is not valid: $infile";
}

if (!$type) {
    die "Type of search not specified. If searching for duplicate elements in same file, please use 'same' as third argument. Use 'other' for search between different files\n";
}

sub help {
  print "

usage:

match_mites_ggsearch_both-ori2.pl  <ggsearch_output_file> <run_name> <E-value_cutoff> <type>

example:
perl match_mites_ggsearch_both-ori2.pl aedes_all_good-tirs-to-self.ggsearch_out aedes 1e-30 same

";
  exit 1;
}

my ($volume, $in_dir, $filename) = File::Spec->splitpath($infile);
my $match_out_path = File::Spec->catpath($volume, $in_dir, $title . "_matches.txt");
my $unique_out_path = File::Spec->catpath($volume, $in_dir, $title . "_unique.txt");
open(my $match_out, ">", $match_out_path);
open(my $unique_out, ">", $unique_out_path);

#print "Opening  search result file . . . ";
my $fa_aln_obj = Bio::SearchIO->new(-format => 'fasta', -file   => $infile);
#print "Opened!\n";
#print Dumper($fa_aln_obj);
#print "\n\n\n";

my %bad_hash;

while( my $result = $fa_aln_obj->next_result ) {
    #print Dumper($result);
    my $query_id = $result->query_name();
    #print "Query:  $query_id\n";
    if ($type eq 'same') {
        if (exists $bad_hash{$query_id}) {
            #print "Query is in bad hash. Continuing to next query.\n";
            next;
        }
        else {
            #print "Seeing query for first time, printing to unique file and adding to bad hash.\n";
            print $unique_out "$query_id\n";
            $bad_hash{$query_id}++;
        }
    }
    while (my $hit = $result->next_hit) {
        my $hsp = $hit->next_hsp;
        my $sig;
        my $hit_id = $hit->name();
        #print "query id: $query_id  hit id: $hit_id\n";
        if (($type eq 'same') and ($query_id eq $hit_id)) {
            #print "Has a hit . . .\t";
            #print "hit is to itself. Skipping.\n";
            next;
        }
        $sig = $hit->significance();
        if (($hsp->percent_identity >= 67.0) and ($sig < $cutoff)) {
            #print "Has a hit . . .\t";
            #print "Hit is good match\n";
            print $match_out "$query_id\t$hit_id\t$sig\n";
            if ($type eq 'same') {
                $bad_hash{$hit_id}++;
            }
        }
        else {
            #print "Has a hit . . .\t";
            #print "Hit is not a good match\n";
            if ($type eq 'other') {
                print $unique_out "$query_id\n";
            }
        }
        if (!$hit and ($type eq 'other')) {
            #print "No hits\n";
            print $unique_out "$query_id\n";
        }
    }
}


close($unique_out);
close($match_out);
