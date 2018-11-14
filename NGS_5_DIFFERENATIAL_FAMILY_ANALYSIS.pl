#!/usr/bin/perl -w
use warnings;
use strict;
use 5.18.0;
use autodie;
use Data::Dumper; $Data::Dumper::Sortkeys = 1;
use Text::Levenshtein::XS qw/distance/;
system("clear");


##### User input #####

my $path_to_file = ""; 

## This is the output from NGS_4_CLEAN_LIBRARIES.pl 
my $file_families = "6_VB_sum_sorted_cleaned.fasta";           


## This are the non-cleaned libraries, that will be clustered and cleaned AFTER comparing to $file_family.
## They come from NGS_2_SUM_SORT.pl output.
my $file_library_one = "9_VB_sum_sorted.fasta";      
my $file_library_two = "12_VB_sum_sorted.fasta";     


## Max. Levenshtein distance to sum up to ancestor. This number should be determined for each experiment by calculating
## the Lv Dist of the Top10/100 to all sequences. Then the histogram will tell you the ideal distance. Usually, the 
## ideal distance lay between 4 and 6. For this analysis the same distance should be used like the on for $file_family.
my $defined_edit_distance = 5; 

## Defines the necessary sum of the readcount in all libraries.
my $readcount_cutoff = 3;


## Defines the nubmer of families that will be displayed, starting with families with biggest changes in either direction.
my $family_cutoff = 5;

## Defines the module (sequence) that is used searhed within the families. Ignore if not used.
my $modul = "CCTACGGGAAAGG";

######################

my (%families, %lib_one, %lib_two);

my ($info, $sequence, $identifier, $readcount, $rpm , $family_index, $clone_number, $counter) = 
   (undef, undef    , undef      , undef     , undef, undef        , undef        , 0);

say "Reading families";
$/= ">";
open (INPUT, "<$path_to_file$file_families");
while (<INPUT>){
    chomp($_);
    next if ($_ eq "");
    if ($_ =~ /(\w+-\d+-\d+.*-\d+\n\w+)/){
        ($info, $sequence)                             = split /\n/, $1;
        ($identifier, $readcount, $rpm, $family_index) = split /-/, $info;
        (undef, $clone_number)                         = split /_K/, $identifier; 

        $families{$counter}{identifier}         = $identifier;
        $families{$counter}{clone_number}       = $clone_number;            
        $families{$counter}{readcount}          = $readcount;
        $families{$counter}{rpm}                = $rpm;
        $families{$counter}{family_index}       = $family_index;
        $families{$counter}{sequence}           = $sequence;
		$families{$counter}{light_module_flag}  = &find_module($sequence);
        $families{$counter}{sum_readcount_lib1} = 0;
        $families{$counter}{sum_readcount_lib2} = 0;
        $families{$counter}{sum_rpm_lib1}       = 0;
        $families{$counter}{sum_rpm_lib2}       = 0;
        $families{$counter}{xfold_lib1}         = "";
        $families{$counter}{xfold_lib2}         = "";
        #$families{$counter}{flag}              = 1;
        last if $counter == 1000;
        $counter++;
    } else {say "ERROR";}
}
close INPUT;
say "File read.";

#print Dumper \%families;


say "Reading library one";
$/= ">";
open (INPUT, "<$path_to_file$file_library_one");
$counter = 0;
while (<INPUT>){
    chomp($_);
    next if ($_ eq "");
    if ($_ =~ /(\w+-\d+-\d+.*\n\w+)/){
        ($info, $sequence)              = split /\n/, $1;
        ($identifier, $readcount, $rpm) = split /-/, $info;
        (undef, $clone_number)          = split /_K/, $identifier; 

        $lib_one{$counter}{identifier}   = $identifier;
        $lib_one{$counter}{clone_number} = $clone_number;            
        $lib_one{$counter}{readcount}    = $readcount;
        $lib_one{$counter}{rpm}          = $rpm;
        $lib_one{$counter}{sequence}     = $sequence;
        $lib_one{$counter}{flag}         = 1;
        last if $counter == 10000;
        $counter++;
    } else {say "ERROR";}
}
close INPUT;

say "Reading library two";
$/= ">";
open (INPUT, "<$path_to_file$file_library_two");
$counter = 0;
while (<INPUT>){
    chomp($_);
    next if ($_ eq "");
    if ($_ =~ /(\w+-\d+-\d+.*\n\w+)/){
        ($info, $sequence)              = split /\n/, $1;
        ($identifier, $readcount, $rpm) = split /-/, $info;
        (undef, $clone_number)          = split /_K/, $identifier; 

        $lib_two{$counter}{identifier}   = $identifier;
        $lib_two{$counter}{clone_number} = $clone_number;            
        $lib_two{$counter}{readcount}    = $readcount;
        $lib_two{$counter}{rpm}          = $rpm;
        $lib_two{$counter}{sequence}     = $sequence;
        $lib_two{$counter}{flag}         = 1;
        last if $counter == 10000;
        $counter++;
    } else {say "ERROR";}
}
close INPUT;

#print Dumper \%lib_two;


say "Calculating Lv Dist for library 1";    
foreach my $key_fam (sort {$a<=>$b} keys %families) {
    foreach my $key_lib (sort {$a<=>$b} keys %lib_one){
        next if !$lib_one{$key_lib}{flag};
        #say $i;
        my $distance = distance($families{$key_fam}{sequence}, $lib_one{$key_lib}{sequence});
        #say $distance if defined($distance);
        if (defined($distance) && $distance <= $defined_edit_distance){        
                #say $current_pool{$i}{identifier}."\t".$distance;
                $lib_one{$key_lib}{flag} = 0;
                $families{$key_fam}{sum_readcount_lib1} += $lib_one{$key_lib}{readcount};
                $families{$key_fam}{sum_rpm_lib1}       += $lib_one{$key_lib}{rpm};
        }
    }
}

say "Calculating Lv Dist for library 2";    
foreach my $key_fam (sort {$a<=>$b} keys %families) {
    foreach my $key_lib (sort {$a<=>$b} keys %lib_two){
        next if !$lib_two{$key_lib}{flag};
        #say $i;
        my $distance = distance($families{$key_fam}{sequence}, $lib_two{$key_lib}{sequence});
        #say $distance if defined($distance);
        if (defined($distance) && $distance <= $defined_edit_distance){        
                #say $current_pool{$i}{identifier}."\t".$distance;
                $lib_two{$key_lib}{flag} = 0;
                $families{$key_fam}{sum_readcount_lib2} += $lib_two{$key_lib}{readcount};
                $families{$key_fam}{sum_rpm_lib2} += $lib_two{$key_lib}{rpm};
        }
    }
}

#print Dumper \%families;

## Maintain %families quality
my $z = 0;
foreach my $key_fam (sort {$a<=>$b} keys %families) {
    if ( ($families{$key_fam}{sum_readcount_lib1} + 
          $families{$key_fam}{sum_readcount_lib2} +
          $families{$key_fam}{readcount}           ) <= $readcount_cutoff ){
        delete $families{$key_fam};
        ++$z;
    }
}
#print Dumper \%families;
say $z;

say "Calculating fold-changes for families";
foreach my $key_fam (sort {$a<=>$b} keys %families) {
    $families{$key_fam}{xfold_lib1} = sprintf "%.2f", $families{$key_fam}{sum_rpm_lib1} / $families{$key_fam}{rpm};
    $families{$key_fam}{xfold_lib2} = sprintf "%.2f", $families{$key_fam}{sum_rpm_lib2} / $families{$key_fam}{rpm};
}

#print Dumper \%families;


foreach my $key_fam (sort {$a<=>$b} keys %families) {

    print $families{$key_fam}{identifier}."\t".$families{$key_fam}{rpm}."\t".$families{$key_fam}{sum_rpm_lib1}."\t".$families{$key_fam}{sum_rpm_lib2}."\t".$families{$key_fam}{light_module_flag}."\t".$families{$key_fam}{sequence}."\n";

}



# geht hoch in 1 (hier einfach die höchsten xfold lib 1 nehmen)
# geht hoch in 2 (hier einfach die höchsten xfold lib 2 nehmen)
# geht runter in 1 (hier muss geschaut werden, dass der wert klein ist und als zweites wird sortiert nach $family{rpm})
# geht runter in 2 (dito)
# geht hoch in 1 && runter in 2 (xfold_lib1 größer 1 und xfold_lib2 kleiner 1, soritert nach xfold_lib1)
# geht hoch in 2 && runter in 1 (dito nur vice versa)



#{
#say "Highest x-fold in lib1";
#my $i = 0;
#my @indices_to_print = sort {$families{$b}{xfold_lib1} <=> $families{$a}{xfold_lib1}           
#                                                       ||
#                             $families{$a}{readcount}  <=> $families{$b}{readcount}}  
#                             keys %families;
#foreach my $index (@indices_to_print) {
#    #say $families{$index}{xfold_lib1};
#    ++$i;
#    last if $i == $family_cutoff;
#}
#}
#
#{
#say "Highest x-fold in lib2";
#my $i = 0;
#my @indices_to_print = sort {$families{$b}{xfold_lib2} <=> $families{$a}{xfold_lib2}           
#                                                       ||
#                             $families{$a}{readcount}  <=> $families{$b}{readcount}}  
#                             keys %families;
#foreach my $index (@indices_to_print) {
#    #say $families{$index}{xfold_lib2};
#    ++$i;
#    last if $i == $family_cutoff;
#}
#}
#
#{
#say "Lowest x-fold in lib1";
#my $i = 0;
#my @indices_to_print = sort {$families{$a}{xfold_lib1} <=> $families{$b}{xfold_lib1}           
#                                                       ||
#                             $families{$b}{sum_readcount_lib2}  <=> $families{$a}{sum_readcount_lib2}}  
#                             keys %families;
#foreach my $index (@indices_to_print) {
#    say $families{$index}{identifier}."\t".
#        $families{$index}{rpm}."\t".
#        $families{$index}{sum_rpm_lib1}."\t".
#        $families{$index}{sum_rpm_lib2}."\t".
#        $families{$index}{xfold_lib1};
#    ++$i;
#    #last if $i == $family_cutoff;
#}
#}


sub find_module {

	my ($seq) = @_;
	use re::engine::TRE max_cost => 2 ;                           
    if ( $seq =~ m/$modul/gi )
	{
		return (1);
	} else
	{
		return (0);
	}     
	no re::engine::TRE;
	
}


say "__END__";

## End of file
## Sorry for the bad code!