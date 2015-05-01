#!/usr/bin/env perl

use strict;
use warnings;
use Bio::AlignIO;
use Bio::SearchIO;
use Bio::Seq;
use Bio::SeqIO;
use Bio::PrimarySeq;
use Bio::Tools::IUPAC;
use FindBin;
use File::Spec;
use Getopt::Long;
use Data::Dumper;

my $verbose;

GetOptions(
  'v' => \$verbose,
);

my $infile = shift @ARGV;

if (!-f $infile) {
  print "Supplied filepath is not valid: $!\n";
  help();
}


sub help {
  print "

usage:

check_query_tirs.pl -v <query_seq_file>

-v verbose output

";
  exit 1;
}

#breakdown directory and filename, create output directory
my ($volume, $in_dir, $filename) = File::Spec->splitpath($infile);


my @fname_fin =  split(/\.fa/, $filename);
my $fname_fin    = $fname_fin[0];
my $out_dir_name = $fname_fin;
my $out_path     = File::Spec->catdir($in_dir, $out_dir_name);

if (!-d $out_path) {
  mkdir($out_path) or die "Can't create dir:$out_path $!\n"; 
}

my $good_out_path = File::Spec->catpath($volume, $out_path, $fname_fin . "_good-tirs.fa");
my $bad_out_path = File::Spec->catpath($volume, $out_path, $fname_fin . "_bad-tirs_trimmed.fa");
my $bad_full_path = File::Spec->catpath($volume, $out_path, $fname_fin . "_bad-tirs_full.fa");

my $good_out = Bio::SeqIO->new(-file => ">$good_out_path", -format => 'Fasta');
my $bad_out = Bio::SeqIO->new(-file => ">$bad_out_path", -format => 'Fasta');
my $bad_full = Bio::SeqIO->new(-file => ">$bad_full_path", -format => 'Fasta');


my $in = Bio::SeqIO->new(-file => $infile, -format => 'Fasta');
my $first_path = File::Spec->catpath($volume, $out_path, "first_half" . ".txt");
my $last_path  = File::Spec->catpath($volume, $out_path, "last_half" . ".txt");

while ( my $seq_obj = $in->next_seq() ) {
    my $seq_id = $seq_obj->id();
    my @short = split( " ", $seq_id);
    my $short_id = $short[0];
    my $seq = $seq_obj->seq();
    my $seq_len = length($seq);
    my $first;
    my $last;
    my $limit = 25;
    my $result = 0;
    my $round = 0;
    
    if ($verbose) {
        print "Seq name: $seq_id\n";
        print "result: $result seq len: $seq_len\n";
    }
    my $seqobj;
    my $alg = "gg";
    $result = search($seq_obj, $limit, $alg, $first_path, $last_path, $good_out, $bad_out, $bad_full, $volume, $out_path, 0);
    if ($verbose) {
        print "first result, gg: $result\n";
    }
    if (!defined $result or $result == 0) {
        $alg = "b2";
        $limit = 175;
        if ($seq_len < 352) {
            $limit = int($seq_len/2);
        }
        $result = search($seq_obj, $limit, $alg, $first_path, $last_path, $good_out, $bad_out, $bad_full, $volume, $out_path, 0);
    }
    elsif ($result == 2 and $round < 2) {
        $round++;
        my $trim_seq = substr($seq, 1);
        $seqobj = Bio::Seq->new(-id => $seq_id, -seq => $trim_seq);
        $result = search($seqobj, $limit, $alg, $first_path, $last_path, $good_out, $bad_out, $bad_full, $volume, $out_path, $result);
        if ($verbose) {
            print "trim2 result: $result\n";
        }
    }
    
    elsif ($result == 3 and $round < 2) {
        $round++;
        my $trim_seq = substr($seq, 0, -1);
        $seqobj = Bio::Seq->new(-id => $seq_id, -seq => $trim_seq);
        $result = search($seqobj, $limit, $alg, $first_path, $last_path, $good_out, $bad_out, $bad_full, $volume, $out_path, $result);
        if ($verbose) {
            print "trim3 result: $result\n";
        }
    }
}



#--------------------------Subroutines---------------------------------#

sub search {
    my $self = shift;    ## seq_obj
    my $limit = shift;
    my $alg = shift;
    my $first_path = shift;
    my $last_path = shift;
    my $good_out = shift;
    my $bad_out = shift;
    my $bad_full = shift;
    my $volume = shift;
    my $out_path = shift;
    my $return;
    
    my $seq_id = $self->id();
    my @short = split( " ", $seq_id);
    my $short_id = $short[0];
    my $seq = $self->seq();
    my $seq_len = length($seq);
    
    my $first = substr($seq, 0, $limit);
    my $last = substr($seq, $seq_len-$limit);
    
    #save the two ends as files to use as inputs for a ggsearch search
    open(my $first_out, ">", $first_path) or die "Error creating $first_path. $!\n";
    open(my $last_out, ">", $last_path) or die "Error creating $last_path. $!\n";
    print $first_out ">first\n$first\n";
    print $last_out ">last\n$last\n";
    close($first_out);
    close($last_out);

    #create fasta output filename then call fasta search program

    my @tir_match_results;
    if ($verbose) {
        print "Alg = $alg\n";
    }
    my $out_opt;
    if ($alg eq "gg") {
        if ($verbose) {
            $out_opt = File::Spec->catpath($volume, $out_path, $short_id . "_ggsearch.out");
        }
        else {
            $out_opt = File::Spec->catpath($volume, $out_path, "ggsearch.out");
        }
        system("ggsearch36 -n -i -T 1 -d 1 $last_path $first_path > $out_opt");
        @tir_match_results = match_tirs($self, $out_opt, 1);
    }
    else {
        if ($verbose) {
            $out_opt = File::Spec->catpath($volume, $out_path, $short_id . "_bl2_seq.out");
        }
        else {
            $out_opt = File::Spec->catpath($volume, $out_path, "bl2_seq.out");
        }
        system("bl2seq -p blastn -W 5 -F F -S 2 -j $last_path -i $first_path -o $out_opt");
        @tir_match_results = match_tirs2($self, $out_opt, 1);
    }
    if ($verbose) {
        print "In search before tir_match Dumper: alg = $alg\n";
        print Dumper(@tir_match_results);
    }
    if ($alg eq "gg") {
        if ($tir_match_results[0][0] == 2) {
            if ($verbose) {
                print "return result:  2\n";
            }
            $return = 2;
        }
        elsif ($tir_match_results[0][0] == 3) {
            if ($verbose) {
                print "return result:  3\n";
            }
            $return = 3;
        }
        else {
            foreach my $row_ref (@tir_match_results) {
                my @entry = @{$row_ref};
                if ($verbose) {
                    print "\nentry[0]: $entry[0]\n";
                    print "entry[1]: $entry[1]\n\n";
                }
                if ($entry[0] == 1) {
                    my %matches = %{$entry[1]};
                    if ($verbose) {
                        print "ggsearch matches Dumper\n";
                        print Dumper(\%matches);
                    }
                    my $pattern = $matches{"hit_seq"} . ".+" . $matches{"query_seq"};
                    if ($verbose) {
                        print "Match pattern = $pattern\n";
                    }
                    $seq =~ m/($pattern)/i;
                    my $trim_seq = $1;
                    my $trim_seq_len = length($trim_seq);
                    
                    if (!defined $trim_seq) {
                        if ($verbose) {
                            print "trim seq not defined\n";
                            print "return result:  0\n";
                        }
                        $return = 0;
                        next;
                    }
                    elsif (abs($seq_len-$trim_seq_len) <= 10) {
                        if ($verbose) {
                            print "printed to good. Seq len: $seq_len  Trim len: $trim_seq_len\n";
                        }
                        if ($seq_len == $trim_seq_len) {
                            $good_out->write_seq($self);
                        }
                        else {
                            my $seqobj = Bio::Seq->new(-id => $seq_id . "_trim", -seq => $trim_seq);
                            $good_out->write_seq($seqobj);
                        }
                        if ($verbose) {
                            print "return result:  1\n";
                        }
                        $return = 1;
                        next;
                    }
                    elsif (abs($seq_len-$trim_seq_len) > 10) {
                        if ($verbose) {
                            print "printed to bad. Seq len: $seq_len  Trim len: $trim_seq_len\n";
                        }
                        $bad_full->write_seq($self);
                        my $seqobj = Bio::Seq->new(-id => $seq_id . "_trim-bad", -seq => $trim_seq);
                        $bad_out->write_seq($seqobj);
                        
                        
                        if ($verbose) {
                            print "return result:  1\n";
                        }
                            $return = 1;
                            next;
                    }
                    
                    else{
                        if ($verbose) {
                            print "return result:  0\n";
                        }
                        $return = 0;
                        next;
                    }
                }
                else {
                    if ($verbose) {
                        print "return result:  0\n";
                    }
                    $return = 0;
                    next;
                }
            }
        }
    }
    else {
        my $num_matches = @tir_match_results;
        
        if (!defined $num_matches) {
            if ($verbose) {
                    print "bl2 match results not defined\n";
                }
            $bad_out->write_seq($self);
            if ($verbose) {
                print "return result:  0\n";
            }
            $return = 0;
        }
        
        elsif ( $num_matches == 0) {
            if ($verbose) {
                print "No bl2 matches.\n";
            }
            $bad_out->write_seq($self);
            if ($verbose) {
                print "return result:  0\n";
            }
            $return = 0;
        }
        
        elsif ($num_matches == 1) {
            my $c;
            foreach my $row_ref (@tir_match_results) {
                $c++;
                my @entry = @{$row_ref};
                if ($verbose) {
                    print "One bl2 match result entry. Current count is $c.  entry[0]:  $entry[0]\n";
                }
                if ($entry[0] == 1) {
                    my %matches = %{$entry[1]};
                    if ($verbose) {
                        print "\n\n\n Entry says match.\n\nbl2seq matches Dumper\n\n\n";
                        print Dumper(\%matches);
                    }
                    my $pattern = $matches{"query_seq"} . ".+" . $matches{"hit_seq"};
                    $seq =~ m/($pattern)/i;
                    my $trim_seq = $1;
                    my $trim_seq_len = length($trim_seq);
                    if ($verbose) {
                        print "Seq length: $seq_len  Trim seq length: $trim_seq_len\n";
                    }
                    
                    if (defined $trim_seq) {
                        if ($trim_seq_len < 50) {
                            if ($verbose) {
                                print "printed to bad\n";
                            }
                            $bad_out->write_seq($self);
                            $return = 1;
                        }
                        elsif (abs($seq_len-$trim_seq_len) <= 10) {
                            if ($verbose) {
                                print "printed to good. Seq len: $seq_len  Trim len: $trim_seq_len\n";
                            }
                            if ($seq_len == $trim_seq_len) {
                                $good_out->write_seq($self);
                            }
                            else {
                                my $seqobj = Bio::Seq->new(-id => $seq_id . "_trim", -seq => $trim_seq);
                                $good_out->write_seq($seqobj);
                            }
                            if ($verbose) {
                                print "return result:  1\n";
                            }
                            $return = 1;
                        }
                        elsif (abs($seq_len-$trim_seq_len) > 10) {
                            if ($verbose) {
                                print "printed to trim. Seq len: $seq_len  Trim len: $trim_seq_len\n";
                            }
                            $bad_full->write_seq($self);
                            my $seqobj = Bio::Seq->new(-id => $seq_id . "_trim-bad", -seq => $trim_seq);
                            $bad_out->write_seq($seqobj);
                            
                            if ($verbose) {
                                print "return result:  1\n";
                            }
                            $return = 1;
                        }

                        else {
                            if ($verbose) {
                                print "printed to bad\n";
                            }
                            $bad_out->write_seq($self);
                            $return = 1;
                        }
                    }
                    else {
                        if ($verbose) {
                            print "printed to bad\n";
                        }
                        $bad_out->write_seq($self);
                        $return = 1;
                    }
                }
            }
        }
        else {
            my $bad = 0;
            my $last_len;
            my $last_trim;
            my $last_eval;
            my $last_hit_pos;
            my $last_query_pos;
            
            
            foreach my $row_ref (@tir_match_results) {
                my @entry = @{$row_ref};
                if ($entry[0] == 1) {
                    my %matches = %{$entry[1]};
                    if ($verbose) {
                        print "bl2seq matches Dumper\n";
                        print Dumper(\%matches);
                    }
                    my $pattern = $matches{"query_seq"} . ".+" . $matches{"hit_seq"};
                    $seq =~ m/($pattern)/i;
                    my $trim_seq = $1;
                    my $trim_seq_len = length($trim_seq);
                    if ($verbose) {
                        print "Seq length: $seq_len  Trim seq length: $trim_seq_len\n";
                    }
                    if (defined $last_len and defined $last_trim) {
                        if ($trim_seq_len == $last_len) {
                            if ($matches{"eval"} < $last_eval) {
                                $last_len = $trim_seq_len;
                                $last_trim = $trim_seq;
                                $last_eval = $matches{"eval"};
                                $last_hit_pos = $matches{"hit"};
                                $last_query_pos = $matches{"query"};
                            }
                        }
                        elsif (abs($trim_seq_len - $last_len) <= 5) {
                            if (abs($matches{"hit"} - $last_hit_pos) < 5 and abs($matches{"query"} - $last_query_pos) < 5) {
                                if ($matches{"eval"} < $last_eval) {
                                    $last_len = $trim_seq_len;
                                    $last_trim = $trim_seq;
                                    $last_eval = $matches{"eval"};
                                    $last_hit_pos = $matches{"hit"};
                                    $last_query_pos = $matches{"query"};
                                }
                            }
                            
                        }
                        elsif ($trim_seq_len > $last_len) {
                            $last_len = $trim_seq_len;
                            $last_trim = $trim_seq;
                            $last_eval = $matches{"eval"};
                            $last_hit_pos = $matches{"hit"};
                            $last_query_pos = $matches{"query"};
                        }
                    }
                    else {
                        $last_len = $trim_seq_len;
                        $last_trim = $trim_seq;
                        $last_eval = $matches{"eval"};
                        $last_hit_pos = $matches{"hit"};
                        $last_query_pos = $matches{"query"};
                        $bad = 1;
                    }
                }
                else {
                    next;
                }
            }
            if ($bad == 0) {
                if ($verbose) {
                    print "printed to bad\n";
                }
                $bad_out->write_seq($self);
                $return = 1;
            }
            else {
                if (defined $last_trim) {
                    if (abs($seq_len-$last_len) <= 10) {
                        if ($verbose) {
                            print "printed to good\n";
                        }
                        if ($seq_len == $last_len) {
                            $good_out->write_seq($self);
                        }
                        else {
                            my $seqobj = Bio::Seq->new(-id => $seq_id . "_trim", -seq => $last_trim);
                            $good_out->write_seq($seqobj);
                        }
                        $return = 1;
                    }
                    else {
                        if ($verbose) {
                            print "printed to bad\n";
                        }
                        $bad_full->write_seq($self);
                        my $seqobj = Bio::Seq->new(-id => $seq_id . "_trim-bad", -seq => $last_trim);
                        $bad_out->write_seq($seqobj);
                        $return = 1;
                    }
                }
                else {
                    if ($verbose) {
                        print "printed to bad\n";
                    }
                    $bad_out->write_seq($self);
                    $return = 1;
                }
            }
        }
    }
    return($return)
}

sub match_tirs {
    my $self       = shift;    ## seq_obj
    my $input_path = shift;
    my $round      = shift;
    
    $self->throw("Need Bio::Seq argument")
        unless ref $self && $self->isa('Bio::Seq');
    if (!-f $input_path) {
        die "Supplied filepath is not valid";
    }
    my $seq = $self->seq();
    $seq =~ s/-//g;
    my $check_seq = $seq;
    $check_seq =~ tr/ATGCatgc/TACGtacg/;
    $check_seq = reverse($check_seq);
    my $seq_len = length($seq);
    my $seq_name = $self->id();
    my $fa_aln_obj = Bio::SearchIO->new(-format => 'fasta', -file => $input_path);
    my @result;
    
    #go through the FASTA input object . . down to the HSP
    while (my $result = $fa_aln_obj->next_result) {
        if ($verbose) {
            print "numHits: ",$result->num_hits,"\n";
        }
        if ($result->num_hits == 0) {
            push @result, [0, [0]];
            last;
        }
        my $hit = $result->next_hit();
        my $hit_name = $hit->name();
        my $hit_len = $result->database_letters();
        while (my $hsp = $hit->next_hsp()) {

            #grab the query, hit, and homology strings
            my $homo_string = $hsp->homology_string;
            my $query_str   = $hsp->query_string;
            my $hit_str     = $hsp->hit_string;
            my $len_homo    = length($homo_string);
            my $len_query   = length($query_str);
            my $len_hit     = length($hit_str);
            #print "homology string\n$homo_string\n";
            my $query_start = $hsp->start('query');
            my $hit_start = $hsp->start('hit');
            if ($verbose) {
                print "query start: $query_start  hit start: $hit_start\nhit string: $hit_str\nquery string: $query_str\n";
            }

            #initialize variables
            my $match_len     = 0;
            my $start_pos     = '';
            my $end_pos       = '';
            my $match_query   = '';
            my $match_hit     = '';
            my $match_mis_aln = 0;
            my $total_mis_aln = 0;
            my $last_good     = 0;
            my $hit_pos;
            my $query_pos;
            #my $match_cutoff = 6;
            my $first_match = 0;
            my $last = 0;
            
            #parse homology string, keeping track of match length and mismatches or gaps
            for (my $count = 0; $count < length($homo_string); $count++) {
                my $homo_char  = substr($homo_string, $count, 1);
                my $query_char = substr($query_str,   $count, 1);
                my $hit_char   = substr($hit_str,     $count, 1);
                if ($round == 1 and $count == 12 and $total_mis_aln >= 7) {
                    if ($match_len <= 3) {
                        $match_len   = 0;
                        $start_pos   = '';
                        $match_query = '';
                        $match_hit   = '';
                        $end_pos     = '';
                        $last = 0;
                        if ($verbose) {
                            print "No TIRs found near start of sequences, resetting counts and ending\n";
                        }
                        last;
                    }
                }
                
                if ($match_len == 0) {
                    #if match length equals 0 and position is not a match, continue to next position
                    if ($homo_char eq " ") {
                        if ($round == 1 or $round == 2) {
                            $total_mis_aln++;
                            next;
                        }
                        else {
                            if ($first_match == 0) {
                                next;
                            }
                            else {
                                $total_mis_aln++;
                                next;
                            }
                        }
                    }
                    
                    #if position is a match, store info, continue to next position
                    elsif ($homo_char eq ":") {
                        $start_pos = $count;
                        $last_good = $count;
                        $match_len++;
                        $match_query .= $query_char;
                        $match_hit   .= $hit_char;
                        $last = 0;
                        if ($first_match == 0) {
                            $first_match = 1;
                        }
                        if ($verbose) {
                            print "Initial match at $start_pos. hit char = $hit_char  \n";
                        }
                        next;
                    }
                }
                if ($match_len == 1) {
                    if ($verbose) {
                        print "match len = 1\n";
                    }
                    if ($homo_char eq " ") {
                        if ($count == 2) {
                            $last_good = $count;
                            $match_len++;
                            $match_query .= $query_char;
                            $match_hit   .= $hit_char;
                            $last = 1;
                            if ($verbose) {
                                print "Mismatch at $count. Allowing. hit char = $hit_char  Length is $match_len\n";
                            }
                            next;
                        }
                        else {
                            $match_len     = 0;
                            $start_pos     = '';
                            $match_query   = '';
                            $match_hit     = '';
                            $match_mis_aln = 0;
                            $last = 0;
                            $total_mis_aln++;
                            if ($verbose) {    
                                print "equal - Mismatch at $count, resetting counts\n";
                            }
                            next;
                        }
                    }
                    else {
                        $last_good = $count;
                        $match_len++;
                        $match_query .= $query_char;
                        $match_hit   .= $hit_char;
                        $last = 0;
                        if ($verbose) {
                            print "Another match at $count. hit char = $hit_char  Length is $match_len\n";
                        }
                        next;
                    }
                }
                elsif ($match_len >= 1 and $match_len < 7) {
                    if ($verbose) {
                        print "Entered match len 1-7 block\n";
                    }
                    
                    #if match length is 1-3 and position is not a match, increment mismatch counter and check if more than one mismatch has occurred
                    if ($homo_char eq " ") {
                        $match_mis_aln++;
                        $total_mis_aln++;
                        
                        if ($count < 6) {
                            if ( $query_char eq "-") {
                                    push @result, [2, [0]];
                                    last;
                                }
                            if ($hit_char eq "-") {
                                push @result, [3, [0]];
                                last;
                            }
                        }
                        
                        #more than one mismatch, reset counters and other info, continue
                        if ($match_len <= 4) {
                            if ($match_mis_aln >= 2 or $last == 1) {
                                $match_len     = 0;
                                $start_pos     = '';
                                $match_query   = '';
                                $match_hit     = '';
                                $match_mis_aln = 0;
                                $last = 0;
                                if ($verbose) {
                                    print "Another Mismatch at $count, resetting counts\n";
                                }
                                next;
                            }
                            else {
                                $match_len++;
                                $last_good = $count;
                                $match_query .= $query_char;
                                $match_hit   .= $hit_char;
                                $last = 1;
                                if ($verbose) {
                                    print "Another Mismatch at $count. hit char = $hit_char  \n";
                                }
                                next;
                            }
                        }
                        elsif ($match_len > 4) {
                            if ($match_mis_aln < 4) {
                                $match_len++;
                                $last_good = $count;
                                $match_query .= $query_char;
                                $match_hit   .= $hit_char;
                                $last = 1;
                                if ($verbose) {
                                    print "Another Mismatch at $count. hit char = $hit_char  \n";
                                }
                                next;
                            }
                            else {
                                $match_len     = 0;
                                $start_pos     = '';
                                $match_query   = '';
                                $match_hit     = '';
                                $match_mis_aln = 0;
                                $last = 0;
                                if ($verbose) {
                                    print "Another Mismatch at $count, resetting counts\n";
                                }
                                next;
                            }
                        }
                        elsif ($total_mis_aln >= 3) {
                            $match_len   = 0;
                            $start_pos   = '';
                            $match_query = '';
                            $match_hit   = '';
                            $end_pos     = '';
                            $last = 0;
                            last;
                        }
                    }
                    
                    #position is a match, store info and continue
                    elsif ($homo_char eq ":") {
                        $last_good = $count;
                        $match_len++;
                        $match_query .= $query_char;
                        $match_hit   .= $hit_char;
                        $last = 0;
                        if ($verbose) {
                            print "Another match at $count. hit char = $hit_char  Length is $match_len\n";
                        }
                        next;
                    }
                }
                
                elsif ($match_len >= 7) {
                    #match length is $match_cutoff or higher. If position is not a match, increment mismatch counter and check if more than 2 mismatches have occurred. If a match, continue.
                    if ($homo_char eq " ") {
                        $match_mis_aln++;
                        $total_mis_aln++;
                        
                        if ($count == length($homo_string)-1) {
                            $end_pos = $last_good;
                            $match_query =~ s/-//g;
                            $match_hit   =~ s/-//g;
                            
                            #reverse complement the match query sequence
                            #$match_query =~ tr/ATGCatgc/TACGtacg/;
                            #$match_query = reverse($match_query);
                            my $match_query_len = length($match_query);
                            my $match_hit_len   = length($match_hit);
                            
                            #find the position in the full sequence of the hit and query match sequences
                            $hit_pos = index(uc($seq), uc($match_hit)) + 1;
                            my $initial_query_pos = index(uc($check_seq), uc($match_query)) + 1;
                            
                            $query_pos = $seq_len - $initial_query_pos;
                            #reverse complement the match query sequence
                            $match_query =~ tr/ATGCatgc/TACGtacg/;
                            $match_query = reverse($match_query);
                            if ($verbose) {
                                print "$seq_name hit_pos:$hit_pos query_pos:$query_pos initial qeury_pos: $initial_query_pos\n";
                            }
                            #store sequence name and the hit and query info
                            my %match = ("hit"   => $hit_pos, "hit_seq" => $match_hit,
                            "query" => $query_pos, "query_seq" => $match_query);
                            
                            push @result, [1, \%match];
                            if ($verbose) {
                                print "Another Mismatch at $count. Reached last loop. Match is long, pushing match info for output\n";
                            }
                            last;
                        }
                        #mismatches under 3, store info and continue
                        if ($match_mis_aln <= 3) {
                            $match_len++;
                            $last_good = $count;
                            $match_query .= $query_char;
                            $match_hit   .= $hit_char;
                            if ($verbose) {
                                print "Another Mismatch at $count, proceeding. hit char = $hit_char\n";
                            }
                            next;
                        }
                        
                        #mismatches 3 or more, store final info for alignment match and end parsing
                        elsif ($match_mis_aln >= 3) {
                            $end_pos = $last_good;
                            $match_query =~ s/-//g;
                            $match_hit   =~ s/-//g;
                            
                            #reverse complement the match query sequence
                            #$match_query =~ tr/ATGCatgc/TACGtacg/;
                            #$match_query = reverse($match_query);
                            my $match_query_len = length($match_query);
                            my $match_hit_len   = length($match_hit);
                            
                            #find the position in the full sequence of the hit and query match sequences
                            $hit_pos = index(uc($seq), uc($match_hit)) + 1;
                            my $initial_query_pos = index(uc($check_seq), uc($match_query)) + 1;
                            
                            $query_pos = $seq_len - $initial_query_pos;
                            #reverse complement the match query sequence
                            $match_query =~ tr/ATGCatgc/TACGtacg/;
                            $match_query = reverse($match_query);
                            if ($verbose) {
                                print "$seq_name hit_pos:$hit_pos query_pos:$query_pos initial qeury_pos: $initial_query_pos\n";
                            }
                            #store sequence name and the hit and query info
                            my %match = ("hit"   => $hit_pos, "hit_seq" => $match_hit,
                            "query" => $query_pos, "query_seq" => $match_query);
                            
                            push @result, [1, \%match];
                            if ($verbose) {
                                print "Another Mismatch at $count. Match is long, pushing match info for output\n";
                            }
                            last;
                        }
                    }
                    
                    #position is a match, store info and continue
                    elsif ($homo_char eq ":") {
                        $last_good = $count;
                        $match_len++;
                        $match_query .= $query_char;
                        $match_hit   .= $hit_char;
                        if ($count == (length($homo_string)-1)) {
                            $end_pos = $last_good;
                            $match_query =~ s/-//g;
                            $match_hit   =~ s/-//g;
                            if ($verbose) {
                                print "Another match at $count. hit char = $hit_char  Length is $match_len\n";
                            }
                            #reverse complement the match query sequence
                            #$match_hit =~ tr/ATGCatgc/TACGtacg/;
                            #$match_hit = reverse($match_hit);
                            my $match_query_len = length($match_query);
                            my $match_hit_len   = length($match_hit);
                            
                            #find the position in the full sequence of the hit and query match sequences
                            $hit_pos = index(uc($seq), uc($match_hit)) + 1;
                            my $initial_query_pos = index(uc($check_seq), uc($match_query)) + 1;
                            $query_pos = $seq_len - $initial_query_pos;
                            if ($verbose) {
                                print "$seq_name hit_pos:$hit_pos query_pos:$query_pos initial_query_pos: $initial_query_pos\n";
                            }
                            
                            $match_query =~ tr/ATGCatgc/TACGtacg/;
                            $match_query = reverse($match_query);
                            #store sequence name and the hit and query info
                            my %match = ("hit"   => $hit_pos, "hit_seq" => $match_hit,
                            "query" => $query_pos, "query_seq" => $match_query);
                            
                            push @result, [1, \%match];
                            last;
                        }
                        if ($verbose) {
                            print "Another match at $count. hit char = $hit_char  Length is $match_len\n";
                        }
                        next;
                    }
                }
            }
            #add in check to see if TIRs were found.
            if (!@result and $end_pos eq '') {
                push @result, [0, [0]];
            }
        }
    }
    if ($verbose) {
        print "result Dumper at end of ggsearch match tirs\n";
        print Dumper(\@result);
    }
    return (@result);
}

sub match_tirs2 {
    my $self       = shift;    ## seq_obj
    my $input_path = shift;
    my $round      = shift;
    
    $self->throw("Need Bio::Seq argument")
        unless ref $self && $self->isa('Bio::Seq');
    if (!-f $input_path) {
        die "Supplied filepath is not valid";
    }
    my $seq = $self->seq();
    $seq =~ s/-//g;
    my $seq_len = length($seq);
    my $seq_name = $self->id();
    my $fa_aln_obj = Bio::SearchIO->new(-format => 'blast', -file => $input_path);
    my @result;
    if ($verbose) {
        print "ggsearch results path: $input_path\n";
    }
    
    #go through the FASTA input object . . down to the HSP
    while (my $result = $fa_aln_obj->next_result) {

        #print "numHits: ",$result->num_hits,"\n";
        if ($result->num_hits == 0) {
            push @result, [0, [0]];
            last;
        }
        my $hit = $result->next_hit();
        my $hit_name = $hit->name();
        my $hit_len = $result->database_letters();
        while (my $hsp = $hit->next_hsp()) {

            #grab the query, hit, and homology strings
            my $homo_string = $hsp->homology_string;
            my $evalue = $hsp->evalue();
            my $query_str   = $hsp->query_string;
            my $hit_str     = $hsp->hit_string;
            my $len_homo    = length($homo_string);
            my $len_query   = length($query_str);
            my $len_hit     = length($hit_str);
            my $query_start = $hsp->start('query');
            my $hit_start = $hsp->start('hit');
            if ($verbose) {
                print "query start: $query_start  hit start: $hit_start\n";
            }

            #initialize variables
            my $match_len     = 0;
            my $start_pos     = '';
            my $end_pos       = '';
            my $match_query   = '';
            my $match_hit     = '';
            my $match_mis_aln = 0;
            my $total_mis_aln = 0;
            my $last_good     = 0;
            my $hit_pos;
            my $query_pos;
            #my $match_cutoff = 6;
            my $first_match = 0;
            
            #parse homology string, keeping track of match length and mismatches or gaps
            for (my $count = 0; $count < length($homo_string); $count++) {
                my $homo_char  = substr($homo_string, $count, 1);
                my $query_char = substr($query_str,   $count, 1);
                my $hit_char   = substr($hit_str,     $count, 1);
                if ($round == 1 and $count == 8 and $total_mis_aln >= 5) {
                    if ($match_len < 3) {
                        $match_len   = 0;
                        $start_pos   = '';
                        $match_query = '';
                        $match_hit   = '';
                        $end_pos     = '';
                        if ($verbose) {
                            print "No TIRs found near start of sequences, resetting counts and ending\n";
                        }
                        last;
                    }
                }
                
                ## skip any seqs that have 2 or more mismatches in the first 3 bases of the TIR
                
                if ($match_len == 0) {
                    if ($count == length($homo_string)-1) {
                        $match_len   = 0;
                        $start_pos   = '';
                        $match_query = '';
                        $match_hit   = '';
                        $end_pos     = '';
                        if ($verbose) {
                            print "Another Mismatch at $count, Reached end. Match is short, ending\n";
                        }
                        last;
                    }
                    
                    #if match length equals 0 and position is not a match, continue to next position
                    elsif ($homo_char eq " ") {
                        if ($round == 1) {
                            $total_mis_aln++;
                            next;
                        }
                        else {
                            if ($first_match == 0) {
                                next;
                            }
                            else {
                                $total_mis_aln++;
                                next;
                            }
                        }
                    }
                    
                    #if position is a match, store info, continue to next position
                    elsif ($homo_char eq "|") {
                        $start_pos = $count;
                        $last_good = $count;
                        $match_len++;
                        $match_query .= $query_char;
                        $match_hit   .= $hit_char;
                        if ($first_match == 0) {
                            $first_match = 1;
                        }
                        if ($verbose) {
                            print "Initial match at $start_pos\n";
                        }
                        next;
                    }
                }
                elsif ($match_len >= 1 and $match_len < 7) {
                    if ($count == length($homo_string)-1 and ($match_len < 6)) {
                        $match_len   = 0;
                        $start_pos   = '';
                        $match_query = '';
                        $match_hit   = '';
                        $end_pos     = '';
                        if ($verbose) {
                            print "Another Mismatch at $count, Reached end. Match is short, ending\n";
                        }
                        last;
                    }
                    
                    #if match length is 1-3 and position is not a match, increment mismatch counter and check if more than one mismatch has occurred
                    elsif ($homo_char eq " ") {
                        $match_mis_aln++;
                        $total_mis_aln++;
                        
                        if ($count == length($homo_string)-1) {
                            $match_len   = 0;
                            $start_pos   = '';
                            $match_query = '';
                            $match_hit   = '';
                            $end_pos     = '';
                            if ($verbose) {
                                print "Another Mismatch at $count, Reached end, $len_homo. Match is short, ending\n";
                            }
                            last;
                        }
                        #allow one mismatch, store info and continue
                        if ($match_mis_aln <= 1) {
                            $match_len++;
                            $last_good = $count;
                            $match_query .= $query_char;
                            $match_hit   .= $hit_char;
                            if ($verbose) {
                                print "First Mismatch at $count\n";
                            }
                            next;
                        }
                        
                        #more than one mismatch, reset counters and other info, continue
                        elsif ($match_mis_aln > 1 and $match_len < 3) {
                            $match_len     = 0;
                            $start_pos     = '';
                            $match_query   = '';
                            $match_hit     = '';
                            $match_mis_aln = 0;
                            if ($verbose) {
                                print "Another Mismatch at $count, resetting counts\n";
                            }
                            next;
                        }
                        elsif ($match_mis_aln < 3 and $match_len >= 3) {
                            $match_len++;
                            $last_good = $count;
                            $match_query .= $query_char;
                            $match_hit   .= $hit_char;
                            if ($verbose) {
                                print "Another Mismatch at $count\n";
                            }
                            next;
                        }
                        elsif ($total_mis_aln >= 3) {
                            $match_len   = 0;
                            $start_pos   = '';
                            $match_query = '';
                            $match_hit   = '';
                            $end_pos     = '';
                            if ($verbose) {
                                print "Another Mismatch at $count, resetting counts and ending\n";
                            }
                            last;
                        }
                    }
                    
                    #position is a match, store info and continue
                    elsif ($homo_char eq "|") {
                        $last_good = $count;
                        $match_len++;
                        $match_query .= $query_char;
                        $match_hit   .= $hit_char;
                        
                        if ($count == length($homo_string)-1) {
                            if ($match_len <= 5) {
                                $match_len   = 0;
                                $start_pos   = '';
                                $match_query = '';
                                $match_hit   = '';
                                $end_pos     = '';
                                if ($verbose) {
                                    print "Another Mismatch at $count, Reached end. Match is short, ending\n";
                                }
                                last;
                            }
                            else {
                                $end_pos = $last_good;
                                $match_query =~ s/-//g;
                                $match_hit   =~ s/-//g;
                                
                                #reverse complement the match query sequence
                                $match_hit =~ tr/ATGCatgc/TACGtacg/;
                                $match_hit = reverse($match_hit);
                                my $match_query_len = length($match_query);
                                my $match_hit_len   = length($match_hit);
                            
                                #find the position in the full sequence of the hit and query match sequences
                                $query_pos = index(uc($seq), uc($match_query)) + 1;
                                $hit_pos = rindex(uc($seq), uc($match_hit)) + 1;
                                
                                $hit_pos = $hit_pos + $match_hit_len;
                                if ($verbose) {
                                    print "$seq_name hit_pos:$hit_pos query_pos:$query_pos $match_hit\n";
                                }
                                #store sequence name and the hit and query info
                                my %match = ("hit"   => $hit_pos, "hit_seq" => $match_hit,
                                "query" => $query_pos, "query_seq" => $match_query, "eval" => $evalue);
                                
                                push @result, [1, \%match];
                                if ($verbose) {
                                    print "Another Mismatch at $count. Match is long, pushing match info for output\n";
                                }
                                last;
                                }
                            }
                        
                        if ($verbose) {
                            print "Another match at $count. Length is $match_len\n";
                        }
                        next;
                    }
                }
                elsif ($match_len >= 7) {
                    
                    #match length is $match_cutoff or higher. If position is not a match, increment mismatch counter and check if more than 2 mismatches have occurred. If a match, continue.
                    if ($homo_char eq " ") {
                        $match_mis_aln++;
                        $total_mis_aln++;
                        
                        #mismatches under 3, store info and continue
                        if ($match_mis_aln <= 3) {
                            $match_len++;
                            $last_good = $count;
                            $match_query .= $query_char;
                            $match_hit   .= $hit_char;
                            if ($verbose) {
                                print "Another Mismatch at $count, proceeding\n";
                            }
                            if ($count == (length($homo_string)-1)) {
                                $end_pos = $last_good;
                                $match_query =~ s/-//g;
                                $match_hit   =~ s/-//g;
                                
                                #reverse complement the match query sequence
                                $match_hit =~ tr/ATGCatgc/TACGtacg/;
                                $match_hit = reverse($match_hit);
                                my $match_query_len = length($match_query);
                                my $match_hit_len   = length($match_hit);
                                
                                #find the position in the full sequence of the hit and query match sequences
                                $query_pos = index(uc($seq), uc($match_query)) + 1;
                                $hit_pos = rindex(uc($seq), uc($match_hit)) + 1;
                                
                                $hit_pos = $hit_pos + $match_hit_len;
                                if ($verbose) {
                                    print "$seq_name hit_pos:$hit_pos query_pos:$query_pos $match_hit\n";
                                }
                                #store sequence name and the hit and query info
                                my %match = ("hit"   => $hit_pos, "hit_seq" => $match_hit,
                                "query" => $query_pos, "query_seq" => $match_query, "eval" => $evalue);
                                
                                push @result, [1, \%match];
                                last;
                            }
                            next;
                        }
                        
                        #mismatches 3 or more, store final info for alignment match and end parsing
                        elsif ($match_mis_aln >= 3) {
                            $end_pos = $last_good;
                            $match_query =~ s/-//g;
                            $match_hit   =~ s/-//g;
                            
                            #reverse complement the match query sequence
                            $match_hit =~ tr/ATGCatgc/TACGtacg/;
                            $match_hit = reverse($match_hit);
                            my $match_query_len = length($match_query);
                            my $match_hit_len   = length($match_hit);
                            
                            #find the position in the full sequence of the hit and query match sequences
                            $query_pos = index(uc($seq), uc($match_query)) + 1;
                            $hit_pos = rindex(uc($seq), uc($match_hit)) + 1;
                            
                            $hit_pos = $hit_pos + $match_hit_len;
                            if ($verbose) {
                                print "$seq_name hit_pos:$hit_pos query_pos:$query_pos $match_hit\n";
                            }
                            #store sequence name and the hit and query info
                            my %match = ("hit"   => $hit_pos, "hit_seq" => $match_hit,
                            "query" => $query_pos, "query_seq" => $match_query, "eval" => $evalue);
                            
                            push @result, [1, \%match];
                            if ($verbose) {
                                print "Another Mismatch at $count. Match is long, pushing match info for output\n";
                            }
                            last;
                        }
                    }
                    
                    #position is a match, store info and continue
                    elsif ($homo_char eq "|") {
                        $last_good = $count;
                        $match_len++;
                        $match_query .= $query_char;
                        $match_hit   .= $hit_char;
                        if ($count == (length($homo_string)-1)) {
                            $end_pos = $last_good;
                            $match_query =~ s/-//g;
                            $match_hit   =~ s/-//g;
                            
                            #reverse complement the match query sequence
                            $match_hit =~ tr/ATGCatgc/TACGtacg/;
                            $match_hit = reverse($match_hit);
                            my $match_query_len = length($match_query);
                            my $match_hit_len   = length($match_hit);
                            
                            #find the position in the full sequence of the hit and query match sequences
                            $query_pos = index(uc($seq), uc($match_query)) + 1;
                            $hit_pos = rindex(uc($seq), uc($match_hit)) + 1;
                            
                            $hit_pos = $hit_pos + $match_hit_len;
                            if ($verbose) {
                                print "$seq_name hit_pos:$hit_pos query_pos:$query_pos $match_hit\n";
                            }
                            #store sequence name and the hit and query info
                            my %match = ("hit"   => $hit_pos, "hit_seq" => $match_hit,
                            "query" => $query_pos, "query_seq" => $match_query, "eval" => $evalue);
                            
                            push @result, [1, \%match];
                            last;
                        }
                        if ($verbose) {
                            print "Another match at $count. Length is $match_len\n";
                        }
                        next;
                    }
                }
            }
            #add in check to see if TIRs were found.
            if ($end_pos eq '') {
                push @result, [0, [0]];
            }
        }
    }
    return (@result);
}
