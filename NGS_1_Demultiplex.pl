#!/usr/bin/perl -w
use warnings;
use strict;
use 5.18.0;
use autodie;
system("clear");


##### User input #####

my $full_path      = "/home/user/";
my $file_name      = "NGS.fastq";

my $five_const     = "GGGCAACTCCAAGCTAGATCTACCGGT";             ## 5' constant region of SELEX pool (5'->3')
my $three_const    = "AAAATGGCTAGCAAAGGAGAAGAACTTTTCACT";	## 3' constant region of SELEX pool (5'->3')

my $length_barcode = 7;

my $barcode_file   = "barcodes.txt";	# Placed in $full_path
my %barcodes       = &get_barcodes($full_path.$barcode_file);

my ($barcodes_recognized, $barcode_correct) = (0, 0);

######################


##### Expected input #####
#
# @title and optional description
# sequence line
# +optional repeat of title line --> in the output here the variable region will be written after the '+'-sign
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


## Reading fastq file

my ($read_counter, $line_counter) = (0, 0,);
my ($identifier, $sequence, $library, $rnd_region, $optional_information, $quality_score);

open (FASTQ_FILE, '<', $full_path.$file_name); 
say "Processing file: $file_name";
while (<FASTQ_FILE>){
    chomp();

    die "Lineissue $_" if ( ($line_counter == 0 ) && (substr($_,0,1) ne "@"));	## First, short quality control

    if      ($line_counter == 0) { $identifier             	= $_; 
    } elsif ($line_counter == 1) { $sequence 			= $_;
    				  ($library, $rnd_region) 	= &qualitycontrol($_);
    } elsif ($line_counter == 2) { $optional_information   	= $_;
    } elsif ($line_counter == 3) { $quality_score  		= $_; 
    } else                       { say "OOPS!";}
    
    ++$line_counter;
    if ( $line_counter > 3) { $line_counter = 0; ++$read_counter; }

    ## Writing output to HD according to barcode
    if ( ($line_counter == 0) && ($rnd_region ne "---") ) {
		my $output = $full_path.$library.".fastq";
   	 	open  (OUTPUT, ">>$output");  
    	print OUTPUT $identifier."\n".$sequence."\n+".$rnd_region."\n".$quality_score."\n";
    	close (OUTPUT);
    }

    #last if $read_counter == 10_000;	## Remove in final run. For debugging only.
    
    my $datestring = localtime();
	say $datestring . ":  " . $read_counter if (($line_counter == 0) && ($read_counter % 100_000 == 0));
}
close FASTQ_FILE;

say "FastQ file read & processed.\n\n";
say "Total number of barcodes recognized:  " . $barcodes_recognized;
say "Total number of correct barcodes:     " . $barcode_correct;
say "Percentage of corrected barcodes:     " . (1-($barcode_correct/$barcodes_recognized))*100 . " %";
say "Percentage sequences recognized:      " . $barcodes_recognized/$read_counter*100 . " %\n\n\n";

#######################



## Subroutines:

## Quality control

sub qualitycontrol {
    my ($seq) = @_;
	use re::engine::TRE (max_cost => 2);	# Allow up to two mismatches while recognizing the constant regions
	if ($seq =~ m/$five_const/){
		return (&assign_library($`), &extract_rnd_region($'));
	} elsif (&rc($seq) =~ m/$five_const/){
		return (&assign_library($`), &extract_rnd_region($'));
	} else {
		# If constant region can not be detected, return "0"
		return ("0", $seq);
	}
	no re::engine::TRE;	
}


## Assigning libraries using error-correctable barcodes

sub assign_library {
    
    my $barcode_to_identify = substr($_[0], -$length_barcode);
	my ($assigned_library, $flag) = (0, 0);

	while (my ($barcode, $library) = each %barcodes) {
		
		# Next line just for statistics
		if ($barcode_to_identify =~ m/$barcode/){
			$barcode_correct++;
		}
		
		use re::engine::TRE (max_cost => 1); 	## Max_cost depends on the designed barcodes
		if ($barcode_to_identify =~ m/$barcode/){
			$assigned_library = $barcodes{$barcode};
			$flag++;
			$barcodes_recognized++;
		}
		no re::engine::TRE;
	}
	return ($assigned_library);
}


## Extract rnd_region

sub extract_rnd_region {
	# For matching the 3' constant region, max_cost is set to a higher value.
	# This increases recovery rate
	use re::engine::TRE (max_cost => 4, max_ins => 1, max_del => 1); 
	if ($_[0] =~ m/$three_const/){
		return ($`);
	} else {
		return ("---");
	}
	no re::engine::TRE;	
}


## Reverse complement

sub rc {
	my $rc = reverse($_[0]);			
	$rc =~ tr/ATGCatgcNn/TACGtacgNn/;	
	return ($rc);
}


## Get barcodes from file and store in hash

sub get_barcodes {
	my %barcodes;
	open (BARCODE_FILE, '<', $_[0]); 
	while (<BARCODE_FILE>){
		# To Do match: "CTAAGTC" => "1"
	}
	return (%barcodes)
}



## End of file
