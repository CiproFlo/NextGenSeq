#!/usr/bin/perl -w
use warnings;
use strict;
use 5.18.0;
use autodie;
#use Data::Dumper;
system("clear");



##### User input #####

my $start_round = 0;
my $end_round   = 37;
my $full_path = "/Volumes/Macintosh\ HD2/__PostDoc\ und\ PhD\ NEU/02_BELEX/_NGS/_seq_data/";

######################


##### Expected input #####
#
# @title and optional description
# sequence line
# +'variable region' ## This one was extracted and then again written on HD by NGS_1_Demultiplex_STREAM_Linux.pl
# quality line
#
# EXAMPLE:
#
# @NS500786:89:HCV7MBGX2:1:11101:24717:1089 1:N:0:ATTACTCG+AGGCTATA
# TATAGTGGATCCGACCGTGGTGCCGTGATCACGGTATCGGATTAGGCCCATACTTATCGCTTTTCTACCTACGTCG
# +
# AAAAAEEEEEEEEAEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEE/EEEEEEEEEE</EEEEEAEAEEEEE/EEE
#
##########################



for (my $round = $start_round; $round <= $end_round; $round++){
    
    ## Some files are not existing, so --> next
    next if ( $round == 24 or $round == 36);

    ## Read input file and sum up during hash generation
    my (%pool, $sequence);             
    my ($read_counter, $line_counter, $totalcount, $diff_seq) = (0, 0, 0, 0);

    my $input_file = $full_path.$round.".fastq";

    open (FASTQ_FILE, '<', $input_file); 
    say "Reading file: $input_file";

    while (<FASTQ_FILE>){
        chomp();

        die "Lineissue $_" if ( ($line_counter == 0 ) && (substr($_,0,1) ne "@"));	## First, short quality control

        if      ($line_counter == 0) {  #$identifier = $_; 
        } elsif ($line_counter == 1) {  #$sequence = $_;
        } elsif ($line_counter == 2) {  $sequence = substr($_, 1, );	 ## this strips the leading '+' sign
                                        
                                        if (!$pool{$sequence}){
                                            $pool{$sequence}{readcount}  = 1;
                                            $pool{$sequence}{rpm}        = 0;
                                            $totalcount++;
                                            $diff_seq++;
                                        } elsif ($pool{$sequence}){
                                            $pool{$sequence}{readcount}++;
                                            $totalcount++;                
                                        }
                                        else {say "ERROR";}

        } elsif ($line_counter == 3) { #$quality_score                     = $_;  # just ignnore
        } else                       { say "OOPS!";}
        
        ++$line_counter;
        if ( $line_counter > 3) { $line_counter = 0; ++$read_counter; }


        #last if $read_counter == 10_000;	## For debugging only.
        my $datestring = localtime();

        say $datestring . ":  " . $read_counter if (($read_counter % 100_000 == 0) && ($line_counter == 0));
    }
    close FASTQ_FILE;

    say "Sequences read and summed up!";
    say "Total number of sequences in library ".$round." : ".$totalcount;
    say "Total number of different sequences : ".$diff_seq."  (".$diff_seq/$totalcount.")";


    ## Calculate RPM foreach sequence
    foreach $sequence (keys %pool){
        $pool{$sequence}{rpm} = sprintf "%.2f", ($pool{$sequence}{readcount}/$totalcount*1_000_000);
    }
    say "RPM calculated!";
    
    
    ## Write output file sorted according to RPM
    my $output = $full_path.$round."_VB_sum_sorted.fasta";
	my $numbering = 1;
    open (OUTPUT, ">$output");
    my @sortsequence = sort { $pool{$b}{rpm} <=> $pool{$a}{rpm} } keys %pool;
    foreach $sequence (@sortsequence) {
        print  OUTPUT ">R".$round."K".$numbering."-".$pool{$sequence}{readcount}."-";
        printf OUTPUT "%.2f", $pool{$sequence}{rpm};
        print  OUTPUT "\n".$sequence."\n";
		$numbering++;
    }
    close OUTPUT;
    say "Output written!";
}

## End of file
