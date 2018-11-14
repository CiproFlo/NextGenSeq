#!/usr/bin/perl -w
use warnings;
use strict;
use 5.18.0;
use autodie;
use Data::Dumper; $Data::Dumper::Sortkeys = 1;
use Statistics::Descriptive;
#use RNA;
use Excel::Writer::XLSX;
system("clear");



##### User input #####

my @rounds_to_process = (37, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13); ## 37 = naive pool!

my $path_to_file    = "/user_folder/seq_data/";               ## Also the directory for the different output files!
my $five_const      = "GGGCAACTCCAAGCTAGATCTACCGGT";          ## Sequence 5' constant region: GGGCAACTCCAAGCTAGATCTACCGGT
my $three_const     = "AAAATGGCTAGCAAAGGAGAAGAACTTTTCACT";    ## Sequence 3' constant region: AAAATGGCTAGCAAAGGAGAAGAACTTTTCACT

my $simple           = 0;                               ## Simple statistic: Total, Unique, Orphan (absolute/percentage)

my $top100           = 0;                               ## Calculates RPM of Top100
my $top_max          = 100;                             ## Top100 can be extended to TopXXX

my $backtrack_top100     = 0;                           ## Backtracking of the Top100 most enriched sequences
my $based_on_library     = 11;                          ## Library from which the Top100 are used for backtracking
my @libraries_to_exclude = (0,);                        ## Libraries NOT inlcuded in analysis e.g. library 0 (= Trash)

my $motifs           = 1;                               ## Searching for predefined motifs
my %motifs = (                                          ## Motifs to search for
    "motif1"   => "AAAAAAAAAAA", 
    "motif2"   => "ATCGATCGATG",
);
my $max_cost   = 0;                                ## Maximum allowd mismatches to the motifs

my @libraries_to_exclude_2 = (0,);               ## Libraries NOT inlcuded in analysis e.g. library 0 (= Trash)


## experimental

my $calculate_mfe      = 0;             ## Calculates MFE with vienna package and gives mean and SD of TopXXX
my $number_mfe         = 100;           ## Take the TopXXX sequences within each pool to calculated MFE statistic
my $standard_parameter = 1;             ## If TRUE, then standard parameter file is used, if note...
my $parameter_file     = "file";        ## specify a valid parameter file, e.g. "path_to_file/rna_andronescu2007.par"

######################



## Main program

my %pools;

foreach my $round (@rounds_to_process){
    
    ## Read input file and put everything into a big hash
    my ($info, $sequence, $identifier, $readcount, $rpm, $clone_number);             
    my $input = $path_to_file.$round."_VB_sum_sorted.fasta";
    say "Reading file: $input";
    $/= ">";
    my $counter = 0;
    open (INPUT, "<$input");
    while (<INPUT>){
        chomp();
        next if ($_ eq "");
        #if ($_ =~ /(\w+-\d+-\d+.*\n\w+)/){
        if (/(\w+-\d+-\d+.*\n\w+)/){
            ($info, $sequence)              = split /\n/, $1;
            ($identifier, $readcount, $rpm) = split /-/, $info;
            next if ($readcount eq "2");

            (undef, $clone_number)          = split /_K/, $identifier; 
            
            $pools{$round}{$counter}{identifier}    = $identifier;
            $pools{$round}{$counter}{clone_number}  = $clone_number;            
            $pools{$round}{$counter}{readcount}     = $readcount;
            $pools{$round}{$counter}{rpm}           = $rpm;
            $pools{$round}{$counter}{sequence}      = $sequence;
            $pools{$round}{$counter}{full_sequence} = $five_const.$sequence.$three_const;
            
            $counter++;
        } else {say "ERROR";}
    }
    close INPUT;
}

say "Hash built.";


## Invocation of subroutines only if user specified it

&simple           if $simple;
&top_hundred      if $top100;
&backtrack        if $backtrack_top100;
&motif_search     if $motifs;
&calculate_mfe    if $calculate_mfe;

say "\nStatistics done!\n\n";


###########################################################

## Subroutines following

sub simple {
    
    ## Task of subroutine "simple":
    ## Count total sequences, unique ones, orphans and sequences at least two times found (total and percentage)
    
    say "Subroutine \"SIMPLE\" is running...";
    my $excel_output = $path_to_file."STATISTIC_SIMPLE.xlsx";
    my $workbook = Excel::Writer::XLSX->new( $excel_output );
    my $worksheet = $workbook->add_worksheet('SIMPLE');
    my $format_header = $workbook->add_format( bold => 1, color => 'white', bg_color => 'grey', align => 'center', text_wrap => 1,);
    my $format_percent = $workbook->add_format(num_format => 0x09 );
    my $col = 0;
    my $row = 0;
    $worksheet->write( $row, 0, "Library"                  , $format_header);
    $worksheet->write( $row, 1, "Total sequences"          , $format_header);
    $worksheet->write( $row, 2, "Unique sequences"         , $format_header);
    $worksheet->write( $row, 3, "Orphans"                  , $format_header);
    $worksheet->write( $row, 4, "Found >1 time (Unique!)"  , $format_header);
    $worksheet->write( $row, 5, "%Orphans"                 , $format_header);
    $worksheet->write( $row, 6, "%Found >1 time"           , $format_header);
    
    foreach my $round (@rounds_to_process){
        my $unique = keys %{ $pools{$round}};
        my ($total_readcount, $orphans) = (0, 0);
        for (my $counter = 0; $counter < keys %{ $pools{$round}}; $counter++){
            $total_readcount += $pools{$round}{$counter}{readcount};
            $orphans++ if $pools{$round}{$counter}{readcount} == 1;
        }
        
        ## Output to EXCEL file
        $row++;
        $worksheet->write( $row, 0, $round);
        $worksheet->write( $row, 1, $total_readcount);
        $worksheet->write( $row, 2, $unique);
        $worksheet->write( $row, 3, $orphans);
        $worksheet->write( $row, 4, ($unique-$orphans));
        $worksheet->write( $row, 5, ($orphans/$total_readcount), $format_percent);
        $worksheet->write( $row, 6, (($total_readcount-$orphans)/$total_readcount), $format_percent);
    }
}


sub top_hundred {

    ## Task of subroutine "Top_hundred":
    ## Sums up RPM of Top100 sequences
    
    say "Subroutine \"TOP_HUNDRED\" is running...";
    
    my $excel_output = $path_to_file."STATISTIC_TOP_HUNDRED.xlsx";
    my $workbook = Excel::Writer::XLSX->new( $excel_output );
    my $worksheet = $workbook->add_worksheet('TOP_HUNDRED');
    my $format_header = $workbook->add_format( bold => 1, color => 'white', bg_color => 'grey', align => 'center', text_wrap => 1,);
    my $format_percent = $workbook->add_format(num_format => 0x09 );
    my $col = 0;
    my $row = 0;
    $worksheet->write( $row, 0, "Library"       , $format_header);
    $worksheet->write( $row, 1, "Sum RPM Top100", $format_header);

    foreach my $round (@rounds_to_process){
        my $sum_rpm_top100;
        for (0..($top_max-1)){ $sum_rpm_top100 += $pools{$round}{$_}{rpm}; }
        $row++;
        $worksheet->write( $row, 0, $round);
        $worksheet->write( $row, 1, $sum_rpm_top100);
    }
}


sub backtrack {

    ## Task of subroutine "backtrack":
    ## Takes Top100 sequences of library $based_on_library and extracts RPM
    ## from everey libray except excluded libraries in @libraries_to_exclude

    say "Subroutine \"BACKTRACK\" is running...";
    my $excel_output = $path_to_file."BACKTRACK.xlsx";
    my $workbook = Excel::Writer::XLSX->new( $excel_output );
    my $worksheet = $workbook->add_worksheet('BACKTRACK');
    my $format_header = $workbook->add_format( bold => 1, color => 'white', bg_color => 'grey', align => 'center', text_wrap => 1,);
    my $format_percent = $workbook->add_format(num_format => 0x09 );
    my $col = 0;
    my $row = 0;
    $worksheet->write( 0, 0, "Library/Identifier" , $format_header);
    my (@sequences_to_track, @sequence_identifier);
    for (0..99){
        $sequences_to_track[$_]  = $pools{$based_on_library}{$_}{sequence};
        $sequence_identifier[$_] = $pools{$based_on_library}{$_}{identifier}
    }
    #print Dumper \@sequences_to_track;
    #print Dumper \@sequence_identifier;    
    for my $sequence_position (0..99){
        ++$col;
        $worksheet->write( $row, $col, $sequence_identifier[$sequence_position], $format_header);
        $worksheet->write( ($#rounds_to_process+1), $col, $sequences_to_track[$sequence_position]);
        foreach my $round (@rounds_to_process){
            next if ( grep( /^$round$/g, @libraries_to_exclude ) );
            $worksheet->write( $round, 0, $round, $format_header);
            for (my $counter = 0; $counter < keys %{ $pools{$round}}; $counter++){
                if ( $sequences_to_track[$sequence_position] eq $pools{$round}{$counter}{sequence} ){
                    $worksheet->write( $round, $col, $pools{$round}{$counter}{rpm});
                    last;
                }
            }
        }    
    }
}


sub motif_search {

    ## Task of subroutine "motif_search":
    ## Search pre-defined motifs or whole sequences (%motifs) within each library
    ## allowing a certain degree of randomization ($max_cost)
    ## Also, if libraries should not be considered, exclude them in @libraries_to_exclude_2
    
    # Was will ich zum Schluss haben: Für jede Runde die aufsummierten RPM und die Namen der Sequenzen mit Sequenz
    
    say "Subroutine \"MOTIF_SEARCH\" is running...";
    my $excel_output = $path_to_file."MOTIF_SEARCH.xlsx";
    my $workbook = Excel::Writer::XLSX->new( $excel_output );
    my $worksheet = $workbook->add_worksheet('MOTIF_SEARCH');
    my $format_header = $workbook->add_format( bold => 1, color => 'white', bg_color => 'grey', align => 'center', text_wrap => 1,);
    my $format_percent = $workbook->add_format(num_format => 0x09 );
    my $col = 0;
    my $row = 1;
    $worksheet->write( 1, 0, "Library/Motif" , $format_header);
    $worksheet->write( 0, 1, "Sum RPM" , $format_header);
    #our @sum_rpm;
    foreach my $key ( sort keys %motifs ){
    	my @sum_rpm;
    	my $row = 1;

        my $current_motif = $motifs{$key};
        ++$col;
        $worksheet->write( 1, $col, $current_motif , $format_header);
        #say $key." / ".$current_motif;
        foreach my $round (@rounds_to_process){
            next if ( grep( /^$round$/g, @libraries_to_exclude_2 ) );
            ++$row;
            $worksheet->write( $row, 0, $round ,);
            #$worksheet->write( $round, 0, $round, $format_header);
            use re::engine::TRE max_cost => $max_cost ;                           
            for (my $counter = 0; $counter < keys %{ $pools{$round}}; $counter++){
                if ( $pools{$round}{$counter}{full_sequence} =~ m/$current_motif/gi ){
                    $sum_rpm[$round] += $pools{$round}{$counter}{rpm};
                }
                #$worksheet->write( $round, $col, $pools{$round}{$counter}{rpm}); 
            }  
            no re::engine::TRE;
            $sum_rpm[$round] =  sprintf "%.2f", $sum_rpm[$round];
            #++$col;
            $worksheet->write( $row, $col, $sum_rpm[$round] ,);

        }
        #@sum_rpm = map {sprintf "%.2f", $_} @sum_rpm;
        #print Dumper \@sum_rpm;
    }
}





## experimental section

sub calculate_mfe {

    ## Task of subroutine "calculate_mfe":
    ## Calculate MFE structure and MFE WITH constant regions
    ## Take MFE, sum up over rounds.
    
    say "Subroutine \"CALCULATE_MFE\" is running...";
    foreach my $round (@rounds_to_process){
        my @mfe;
        my $counter; 
        my @clone_numbers = sort {$a<=>$b} keys %{ $pools{$round} };  ## falsch!!!!!!! 
        if (!$standard_parameter){RNA::read_parameter_file($parameter_file);}
        foreach my $clone_number (@clone_numbers){
            my $sequence = $five_const.$pools{$round}{$clone_number}{sequence}.$three_const;
            (undef, my $mfe) = RNA::fold($sequence);
            push (@mfe, $mfe);
            ++$counter;
            last if $counter == $number_mfe;
            }
            my $stat = Statistics::Descriptive::Full->new();
            $stat->add_data(@mfe) ;
            say "Round: ".$round;
            print  "Number calculated MFEs: "        . $stat->count()        ."\n";
            printf "Mean MFE    %.2f\t(kcal/mol)\n"  , $stat->mean()              ;
            printf "SD   MFE    ±%.2f\t(kcal/mol)\n" , $stat->standard_deviation();
            printf "Min  MFE    %.2f\t(kcal/mol)\n"  , $stat->min()               ;
            printf "Max  MFE    %.2f\t(kcal/mol)\n"  , $stat->max()               ;
            print  "-"x20;
    }
}

## End of file
