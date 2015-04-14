#!/usr/bin/perl

use strict;
use warnings;

use Bio::SearchIO;
use File::Spec;
use Data::Dumper;
use FindBin;

#make sure only one input file is specified
my $len = @ARGV;
if ($len != 2) {
    die "Please specify an input file and a run name.\n";
}
my $infile = shift @ARGV;
my $title = shift @ARGV;

if (!-f $infile) {
    die "Supplied filepath is not valid: $infile";
}

my %family_to_super = ("Mirage" => "CMC", "CACTA" => "CMC", "Chapaev" => "CMC", "ISL2EU" => "PIF_harbinger", "Tc1-Mariner" => "Tc_mariner", "MuDR" => "MULE", "Ginger2_TDD" => "Ginger2", "PIF-Harbinger" => "PIF_harbinger", "Ginger2/TDD" => "Ginger2");

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
my $range_out_path = File::Spec->catpath($volume, $in_dir, $title . "_range.txt");
my $long_partial_path = File::Spec->catpath($volume, $in_dir, $title . "_long.txt");
my $short_patial_path = File::Spec->catpath($volume, $in_dir, $title . "_short.txt");
my $no_matches_out_path = File::Spec->catpath($volume, $in_dir, $title . "_no-matches.txt");
my $retros_out_path = File::Spec->catpath($volume, $in_dir, $title . "_retros.txt");
my $heli_out_path = File::Spec->catpath($volume, $in_dir, $title . "_heli-mav.txt");

open(my $full_out, ">", $full_out_path);
open(my $range_out, ">", $range_out_path);
open(my $long_out, ">", $long_partial_path);
open(my $short_out, ">", $short_patial_path);
open(my $no_matches_out, ">", $no_matches_out_path);
open(my $retros_out, ">", $retros_out_path);
open(my $heli_out, ">", $heli_out_path);


print "Opening  search result file . . . ";
my $fa_aln_obj = Bio::SearchIO->new(-format => 'fasta', -file   => $infile);
print "Opened!\n";
my %seen_hash;
my $c;

RESULT: while( my $result = $fa_aln_obj->next_result ) {
    my %results_hash;
    #print Dumper($result);
    my $query_id = $result->query_name();
    if (exists $seen_hash{$query_id}) {
        next RESULT;
    }
    $c++;
    #print "Processing query sequence $c\n";
    $seen_hash{$query_id}++;
    my $query_len = $result->query_length();
    HIT: while ( my $hit = $result->next_hit ) {
        my %track_hash;
        my @family;
        my $fam;
        my $hsp_count = 0;
        HSP: while (my $hsp = $hit->next_hsp ) {
            #print "Has a hit . . .\t";
            my $hit_id = $hit->name();
            my $query_strand = $hsp->strand('query');
            my $hit_strand = $hsp->strand('hit');
            if ($hsp_count > 0) {
                if ($query_strand ne $track_hash{"query_strand"} or $hit_strand ne $track_hash{"hit_strand"}) {
                    next HSP;
                }
                else {
                    $hsp_count++;
                }
                
            }
            else {
                $track_hash{"query_strand"} = $query_strand;
                $track_hash{"hit_strand"} = $hit_strand;
                $hsp_count++;
            }
            if ($hit_id =~ m/ClassI:/) {
                print $retros_out "$query_id\t$hit_id\n";
                next RESULT;
            }
            elsif ($hit_id =~ m/Maverick|Helitron/i) {
                print $heli_out "$query_id\t$hit_id\n";
                next RESULT;
            }
            elsif ($hit_id =~ m/Banshee/) {
                $fam = "Banshee";
            }
            else {
                @family = split(":", $hit_id);
                $fam = $family[-1];
            }
            if ($fam eq "Sola") {
                #print "Hit is to Sola . . .";
                if ($hit_id =~ m/Sola1/) {
                    $fam = "Sola1";
                    #print "Sola1.\n";
                }
                elsif ($hit_id =~ m/Sola-1/) {
                    $fam = "Sola1";
                    #print "Sola1.\n";
                }
                elsif ($hit_id =~ m/Sola2/) {
                    $fam = "Sola2";
                    #print "Sola2.\n";
                }
                elsif ($hit_id =~ m/Sola3/) {
                    $fam = "Sola3";
                    #print "Sola3.\n";
                }
                else {
                    print "Error, Sola element found but could not determine which type.\n";
                }
            }
            
            my $super = $fam;
            if (exists $family_to_super{$fam}) {
                $super = $family_to_super{$fam};
            }
            my $hit_len = $hit->length;
            my $sig = $hit->significance();
            my $percentid = $hsp->percent_identity();
            my $hit_start = $hsp->start('hit');
            my $hit_end = $hsp->end('hit');
            my $match_len = ($hit_end - $hit_start) + 1;
            my $len_ratio = int((($match_len/$hit_len)*100)+.5);
            
            
            if ($hsp_count > 1) {
                $len_ratio = $len_ratio + $track_hash{"len"};
            }
            $track_hash{"info"} = "$query_id\t$query_len\t$hit_id\t$super\t$len_ratio\t$hit_len\t$hit_start-$hit_end\t$sig\t$percentid\n";
            $track_hash{"len"} = $len_ratio;
            $track_hash{"super"} = $super;
            next;
        }
        #print Dumper(\%element_char_hash);
        my $super = $track_hash{"super"};
        print "query len: $query_len  Super: $super\n";
        if ($track_hash{"len"} >= 85) {
            if (exists $results_hash{"full"}) {
                if ($results_hash{"full"}{"len"} < $track_hash{"len"}) {
                    $results_hash{"full"}{"info"} = $track_hash{"info"};
                    $results_hash{"full"}{"len"} = $track_hash{"len"};
                    next HIT;
                }
                else {
                    next HIT;
                }
            }
            else {
                $results_hash{"full"}{"info"} = $track_hash{"info"};
                $results_hash{"full"}{"len"} = $track_hash{"len"};
                next HIT;
            }
        }
        elsif ($track_hash{"len"} >= 50 and $track_hash{"len"} < 85) {
            if (exists $results_hash{"full"}) {
                next HIT;
            }
            elsif (exists $results_hash{"range"}) {
                if ($track_hash{"len"} > $results_hash{"range"}{"len"}) {
                    $results_hash{"range"}{"info"} = $track_hash{"info"};
                    $results_hash{"range"}{"len"} = $track_hash{"len"};
                    next HIT;
                }
                else {
                    next HIT;
                }
            }
                    
            elsif ($super ne "?") {
                if ($query_len >= $element_char_hash{$super}{'min'} && $query_len <= $element_char_hash{$super}{'max'}) {
                    $results_hash{"range"}{"info"} = $track_hash{"info"};
                    $results_hash{"range"}{"len"} = $track_hash{"len"};
                    next HIT;
                }
                else {
                    $results_hash{"long"}{"info"} = $track_hash{"info"};
                    $results_hash{"long"}{"len"} = $track_hash{"len"};
                    next HIT;
                }
            }
            else {
                $results_hash{"long"}{"info"} = $track_hash{"info"};
                $results_hash{"long"}{"len"} = $track_hash{"len"};
                next HIT;
            }
        }
        else {
            if (exists $results_hash{"short"}) {
                if ($track_hash{"len"} > $results_hash{"short"}{"len"}) {
                    $results_hash{"short"}{"info"} = $track_hash{"info"};
                    $results_hash{"short"}{"len"} = $track_hash{"len"};
                    next HIT;
                }
                else {
                    next HIT;
                }
            }
            else {
                $results_hash{"short"}{"info"} = $track_hash{"info"};
                $results_hash{"short"}{"len"} = $track_hash{"len"};
                next HIT;
            }
        }
    }
    if (exists $results_hash{"full"}) {
        print $full_out "$results_hash{'full'}{'info'}";
    }
    elsif (exists $results_hash{"range"}) {
        print $range_out "$results_hash{'range'}{'info'}";
    }
    elsif (exists $results_hash{"long"}) {
        print $long_out "$results_hash{'long'}{'info'}";
    }
    elsif (exists $results_hash{"short"}) {
        print $short_out "$results_hash{'short'}{'info'}";
    }
    else {
        print $no_matches_out "$query_id\n";
        next;
    }
}

close($full_out);
close($range_out);
close($long_out);
close($short_out);
close($no_matches_out);
close($retros_out);
close($heli_out);
