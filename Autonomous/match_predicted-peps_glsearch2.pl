#!/usr/bin/perl

use strict;
use warnings;

use Bio::SearchIO;
use File::Spec;
use Data::Dumper;
use FindBin;

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

my %family_to_super = {"Mirage" => "CMC", "CACTA" => "CMC", "Chapaev" => "CMC", "ISL2EU" => "PIF/harbinger", "Tc1-Mariner" => "Tc/mariner", "MuDR" => "MULE", "Ginger2/TDD" => "Ginger2"};

print "\nProcessing TE TPase length table\n";
my $TE_char_path = $FindBin::Bin . "/DNA-TE_TPase-len.table";
open(my $TE_in, "<", $TE_char_path) or die "Error reading $TE_char_path . $!";
my %element_char_hash;

while (my $line = <$TE_in>) {
  chomp $line;
  my @split = split("\t", $line);
  my $ele_name = $split[0];
  $element_char_hash{$ele_name}{'min'} = $split[1];
  $element_char_hash{$ele_name}{'max'} = $split[2];
}
close($TE_in);


my ($volume, $in_dir, $filename) = File::Spec->splitpath($infile);
my $full_out_path = File::Spec->catpath($volume, $in_dir, $title . "_full.txt");
my $range_out_path = File::Spec->catpath($volume, $in_dir, $title . "_in_range.txt");
my $long_partial_path = File::Spec->catpath($volume, $in_dir, $title . "_long-partial.txt");
my $short_patial_path = File::Spec->catpath($volume, $in_dir, $title . "_short-partial.txt");
my $no_matches_out_path = File::Spec->catpath($volume, $in_dir, $title . "_no-matches.txt");
my $retros_out_path = File::Spec->catpath($volume, $in_dir, $title . "_retros.txt");

open(my $full_out, ">", $full_out_path);
open(my $range_out, ">", $range_out_path);
open(my $long_partial_out, ">", $long_partial_path);
open(my $short_partial_out, ">", $short_patial_path);
open(my $no_matches_out, ">", $no_matches_out_path);
open(my $retros_out, ">", $retros_out_path);


print "Opening  search result file . . . ";
my $fa_aln_obj = Bio::SearchIO->new(-format => 'fasta', -file   => $infile);
print "Opened!\n";

RESULT: while( my $result = $fa_aln_obj->next_result ) {
    my %results_hash;
    #print Dumper($result);
    my $query_id = $result->query_name();
    my $query_len = $result->query_length();
    #print "Query:  $query_id\t$query_len\n";
    HIT: while ( my $hit = $result->next_hit ) {
        my %track_hash;
        my $hsp_count = 0;
        while (my $hsp = $hit->next_hsp ) {
            $hsp_count++;
            #print "Has a hit . . .\t";
            my $hit_id = $hit->name();
            if ($hit_id =~ m/ClassI:/) {
                print $retros_out "$query_id\t$hit_id\n";
                next RESULT;
            }
            my @fam = split(":", $hit_id);
            my $fam = $fam[-1];
            my $super = $fam;
            if (exists $family_to_super{$fam}) {
                $super = $family_to_super{$fam};
            }
            my $hit_len = $hit->length;
            my $sig = $hit->significance();
            my $hit_start = $hsp->hit_start();
            my $hit_end = $hsp->hit_end();
            my $match_len = $hit_end - $hit_start + 1;
            my $len_ratio = int((($hit_len/$match_len)*100) + .5);
            if ($hsp_count > 1) {
                $len_ratio = $len_ratio + $track_hash{"len"};
            }
            if ($len_ratio >= 85) {
                if (exists $track_hash{"info"}) {
                    if ($track_hash{"len"} < $len_ratio) {
                        $track_hash{"info"} = "$query_id\t$hit_id\t$len_ratio\n";
                        $track_hash{"len"} = $len_ratio;
                        next;
                    }
                    else {
                        next;
                    }
                }
                else {
                    #print "Hit is good match and possibly full length\n";
                    $track_hash{"info"} = "$query_id\t$hit_id\t$len_ratio\n";
                    $track_hash{"len"} = $len_ratio;
                    next;
                }
            }
            elsif ($len_ratio >= 50 and $len_ratio < 85) {
                if (exists $track_hash{"info"}) {
                    if ($query_len >= $element_char_hash{$super}{'min'} && $query_len <= $element_char_hash{$super}{'max'}) {
                        $track_hash{"info"} = "$query_id\t$hit_id\t$len_ratio\n";
                        $track_hash{"len"} = $len_ratio;
                        next;
                    }
                    else {
                        $track_hash{"info"} = "$query_id\t$hit_id\t$len_ratio\n";
                        $track_hash{"len"} = $len_ratio;
                        next;
                    }
                }
                else {
                    if ($query_len >= $element_char_hash{$super}{'min'} && $query_len <= $element_char_hash{$super}{'max'}) {
                        $track_hash{"info"} = "$query_id\t$hit_id\t$len_ratio\n";
                        $track_hash{"len"} = $len_ratio;
                        next;
                    }
                    else {
                        $track_hash{"info"} = "$query_id\t$hit_id\t$len_ratio\n";
                        $track_hash{"len"} = $len_ratio;
                        next;
                    }
                }
            }
            else {
                #print "Hit is good match but very short\n";
                $track_hash{"info"} = "$query_id\t$hit_id\t$len_ratio\n";
                $track_hash{"len"} = $len_ratio;
                next;
            }
        }
    }
}

close($full_out);
close($range_out);
close($long_partial_out);
close($short_partial_out);
close($no_matches_out);
