#!/usr/bin/perl -w
use warnings;
use strict;
use 5.18.0;
use autodie;
use Data::Dumper; 
$Data::Dumper::Sortkeys = 1;
use Text::Levenshtein::XS qw/distance/;
system("clear");


##### User input #####

my $path_to_file          = "";       ## Also the directory for the different output files!

my @libraries_to_clean    = (6, 9, 12,);                                  ## 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13);

my $defined_edit_distance = 5;                                     ## Max. Levenshtein distance to sum up to ancestor. This number  
                                                                   ## should be determined for each experiment by calculating
                                                                   ## the Lv Dist of the Top10/100 to all sequences. Then the
                                                                   ## histogram will tell you the ideal distance. Usually, the 
                                                                   ## ideal distance lay between 4 and 6.                  
######################


foreach my $round (@libraries_to_clean) {
    say "Processing round $round";
    my %current_pool;
    my ($info, $sequence, $identifier, $readcount, $rpm, $clone_number, $counter) = ( undef, undef, undef, undef, undef, undef, 0);
    my $input = $path_to_file.$round."_VB_sum_sorted.fasta";
    say "Reading file: $input";
    $/= ">";
    open (INPUT, "<$input");
    while (<INPUT>){
        chomp($_);
        next if ($_ eq "");
        if ($_ =~ /(\w+-\d+-\d+.*\n\w+)/){
            ($info, $sequence)              = split /\n/, $1;
            ($identifier, $readcount, $rpm) = split /-/, $info;
            (undef, $clone_number)          = split /_K/, $identifier; 
            
            $current_pool{$counter}{identifier}   = $identifier;
            $current_pool{$counter}{clone_number} = $clone_number;            
            $current_pool{$counter}{readcount}    = $readcount;
            $current_pool{$counter}{rpm}          = $rpm;
            $current_pool{$counter}{sequence}     = $sequence;
            $current_pool{$counter}{flag}         = 1;

            $counter++;
        } else {say "ERROR";}
    }
    close INPUT;
    say "File read.";
    
    #print Dumper \%current_pool;
    
    my %cleaned_pool;
    
    foreach my $key (sort {$a<=>$b} keys %current_pool) {
        next if !$current_pool{$key}{flag};
        ## Now copy the ancestor into the cleaned HoH
        $cleaned_pool{$key}{identifier}   = $current_pool{$key}{identifier};
        $cleaned_pool{$key}{clone_number} = $current_pool{$key}{clone_number};            
        $cleaned_pool{$key}{readcount}    = $current_pool{$key}{readcount};
        $cleaned_pool{$key}{rpm}          = $current_pool{$key}{rpm};
        $cleaned_pool{$key}{sequence}     = $current_pool{$key}{sequence};
        say $cleaned_pool{$key}{identifier}." in cleaned";
        for (my $i = ($key + 1); $i < $counter; ++$i){
            #say $i;
            my $distance = distance($current_pool{$key}{sequence}, $current_pool{$i}{sequence}) if $current_pool{$i}{flag};
            #say $distance if defined($distance);
            if (defined($distance) && $distance <= $defined_edit_distance){        
                    say $current_pool{$i}{identifier}."\t".$distance;
                    $current_pool{$i}{flag} = 0;
                    $cleaned_pool{$key}{readcount} += $current_pool{$i}{readcount};
                    $cleaned_pool{$key}{rpm}       += $current_pool{$i}{rpm};
            }
        }
    }

    my $output = $path_to_file.$round."_VB_sum_sorted_cleaned.fasta";
    my $family_index = 1;
    open (OUTPUT, ">$output");
    my @indices = sort {$cleaned_pool{$b}{rpm}          <=>  $cleaned_pool{$a}{rpm}            
                                                 ||
                        $cleaned_pool{$a}{clone_number}  <=>  $cleaned_pool{$b}{clone_number}}  
                        keys %cleaned_pool;
    foreach my $counter (@indices) {
        print  OUTPUT ">".$cleaned_pool{$counter}{identifier}."-".$cleaned_pool{$counter}{readcount}."-";
        printf OUTPUT "%.2f", $cleaned_pool{$counter}{rpm};
        print  OUTPUT "-".$family_index."\n".$cleaned_pool{$counter}{sequence}."\n";
        ++$family_index;
    }
    close OUTPUT;
    say "Output written!";


}

say "__END__";

## End of file