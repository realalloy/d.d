#!/usr/bin/perl -w
#*******************************************************************************
#                                 profil.pl                                    *
#                 Copyright (C) Nicolas Le Novère and Marine Dumousseau 2009   *
# Run the program melting iteratively on a nucleic acid sequence entered from  *
#  stdin (it can be a file redirected with 'multi.pl < file.seq')              *
#*******************************************************************************/
#
#      This program is free software; you can redistribute it and/or modify
#      it under the terms of the GNU General Public License as published by
#      the Free Software Foundation; either version 2 of the License, or
#      (at your option) any later version.
#
#      This program is distributed in the hope that it will be useful,
#      but WITHOUT ANY WARRANTY; without even the implied warranty of
#      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#      GNU General Public License for more details.
#
#      You should have received a copy of the GNU General Public License
#      along with this program; if not, write to the Free Software
#      Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
#
#      Nicolas Le Novère and Marine Dumousseau
#      EMBL-EBI, Wellcome-Trust Genome Campus
#      Hinxton Cambridge, CB10 1SD, UK
#      lenov@ebi.ac.uk

       ###################################################################
       # Usage is: profil.pl -Iinfile -Wwindow < inputfile > outputfile  #
       # Where inputfile contains one sequence per line. No space before #
       # sequence and at least one space before extra information        # 
       # Write the parameters of melting in an input file. See manual.   #
       ###################################################################

use strict;

my $VERSION = 2;

my $argument;          # one of the arguments
my $infile = "infile"; # contains the parameters of the run except the sequence
my $window = 10;       # contains the length of the window to analise
my %nucleic_acid ;     # the nucleic acid of the analysis

##########################
# Processes the arguments 
##########################

if (not @ARGV){
    usage();
    exit();
}

if (join("",@ARGV) =~ /-H/i){
    usage();
    exit();
}
foreach $argument (@ARGV){
    if ($argument =~ /^-I/i){        # infile specification
	$argument =~ /^-I([\w\.]*)/i;
	if (defined $1){$infile = $1;}
	else {warning_infile();}
    } elsif ($argument =~ /^-W/i){   # window width specification
	$argument =~ /^-W(\d*)/i;
	if (defined $1){$window = $1;}
	else {warning_window();}
    } else {
	print "Oups! I did not recognise the option $argument.\n";
	print "I do not take it into account\n";
    }
}

#################################
# Reads the sequences to analyse
#################################
# it would be desirable to treat a multisequence FASTA file as
# a set of independant sequences rather than to concatenate 
# everything.


while (<STDIN>){
    if (/^\s*>/){next;}            # remove the fasta info line
    chomp();                       # remove end of line
    s/#.*//;                       # remove the comments
    s/[^AGCTU]//g;                 # keeps only AGCTU, case insensitive
    $_ = uc($_);                   # switch everything uppercase
    $nucleic_acid{"sequence"}.=$_; # append
}
#print "DEBUG--> sequence is: ",$nucleic_acid{"sequence"},"\n";
# Note that each line could contain other elements used in derived programs, but separed by spaces

#############
# Here we go
#############

@{$nucleic_acid{"results"}}=compute_tm($nucleic_acid{"sequence"},$window);
#print "DEBUG--> results: ",$nucleic_acid{"results"},"\n";
print_result();

sub usage{
    print <<EOU;
Usage is: profil.pl [-H] -Iinfile -Wwindow < inputfile > outputfile
where   
    -H, -h, -help print this message;
    -Iinfile      file containing the parameters
    -Wwindow      width of the window analysed at each run
    inputfile     file containing the sequence
    outputfile    file containing the results
EOU
}

sub warning_infile{
    print <<EOWI;
Except the sequences, all parameters have to be contained a configuration
       file, and the script run as:
       prompt> ./profil.pl -Iconfig_file -Wwindow < inputfile > outputfile)
Since no configuration file has been provided (option -I), the file infile
is assumed. See the user-guide of melting for the format of this file.
EOWI
}

sub warning_window{
    print "Since no window length specification has been entered (option -W), a default\n";
    print "length of 10 nucleotides is assumed\n";
}

sub compute_tm{
    my ($sequence,$window)=@_;
    my $i;                     # loop counter
    my @results;               # array of hashes containing the results of one nucleic acid
    
    my $seqlength = length $sequence;
#    print "DEBUG--> length of the sequence is: ",$seqlength,"\n";
    if ($seqlength < $window){$window = $seqlength;} # in case of very short sequences
    
    my $half_window = sprintf ("%d",$window / 2);    # absolute value. 
#    print "DEBUG--> half of the window is: ",$half_window,"\n";
    # Tm is reported to the middle of the string, we have to populate the first half
    # Note the '<', index beginning at 0
    for ($i=0 ; $i<$half_window ; $i++){
	$results[$i]->{"subsequence"} = 'X'x$window;   
	$results[$i]->{"enthalpy"}    = "000000";   
       	$results[$i]->{"entropy"}     = "000.00";   
	$results[$i]->{"tm"}          = "-274";   
#	print "DEBUG--> results[",$i,"] is: ",%{$results[$i]},"\n";
    }
    
    # Now we begin the actual analysis. The sequence move from the beginning, but the
    # temperatures from the middle of the first substring
    my $offset = 0;
    while (($offset+$window) <= $seqlength){
	my $subsequence = substr($sequence,$offset,$window);
#	print "DEBUG--> subsequence is: ",$subsequence,"\n";
	$results[$half_window+$offset]->{"subsequence"} = $subsequence;   
	my @rawresults = `melting -I$infile -S$subsequence -q`;
	foreach (@rawresults){
	    if (/Enthalpy/ ){
		($results[$half_window+$offset]->{"enthalpy"}) = (split)[1];   
	    } elsif (/Entropy/ ){
		($results[$half_window+$offset]->{"entropy"})  = (split)[1];   
	    } elsif (/Melting/ ){
		($results[$half_window+$offset]->{"tm"})       = (split)[2];   
	    } else { # do nothing, this is not suppose to occur but ...
	    } 
	}
#	print "DEBUG--> results[",$half_window+$offset,"] is: ",%{$results[$half_window+$offset]},"\n";
	$offset++;
    }
	
    # fill the remaining position of the array
    # Note that $offset is incremented just before to quit the loop,
    # hence the $half_window+$offset, identical to that of the while loop.
    
    for ( $i = $half_window + $offset ; $i < $seqlength ; $i++ ){ 
	$results[$i]->{"subsequence"} = 'X'x$window;   
	$results[$i]->{"enthalpy"}    = "000000";   
       	$results[$i]->{"entropy"}     = "000.00";   
	$results[$i]->{"tm"}          = "-274";   
#	print "DEBUG--> results[",$i,"] is: ",%{$results[$i]},"\n";
    }
    return @results;
}

sub print_result{
    my $i;         # loop counter
    print ("Base,Enthalpy,Entropy,Tm\n");
    my $seqlength = length $nucleic_acid{"sequence"};
#    print "DEBUG--> length of the sequence is: ",$seqlength,"\n";
    for ($i=0 ; $i<$seqlength ; $i++){
	printf "%s\t", $nucleic_acid{"results"}[$i]->{"subsequence"};
	printf "%7d\t    ",$nucleic_acid{"results"}[$i]->{"enthalpy"};
	printf "%7.2f\t",$nucleic_acid{"results"}[$i]->{"entropy"};
	printf "%7.2f    ",$nucleic_acid{"results"}[$i]->{"tm"};
	print "\n";
    }
}

