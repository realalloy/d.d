eval '(exit $?0)' && eval 'exec perl -S $0 ${1+ "$@"}' && eval 'exec perl -S $0 $argv:q' if 0;

##############################################################################
#                                   MELTING                                  #
# This program   computes for a nucleotide probe, the enthalpie, the entropy #
# of the helix-coil transition, and then its melting temperature.            #
# Three types of hybridisation are possible: DNA/DNA, DNA/RNA, and RNA/RNA.  #
#         Copyright (C) Nicolas Le Novère and Marine Dumousseau 1997-2013    #
##############################################################################
#     This program is free software; you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation; either version 2 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program; if not, write to the Free Software
#     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA#
#
#     Nicolas Le Novere
#     Babraham Institute, Babraham Research Campus
#     Babraham CB22 3AT Cambridge United-Kingdom.
#     n.lenovere@gmail.com
#      
#     Marine Dumousseau
#     EMBL-EBI, Wellcome-Trust Genome Campus
#     Hinxton CB10 1SD Cambridge United-Kingdom. 
#     marine@ebi.ac.uk  


use strict;
use Tk;
require Tk::LabEntry;
require Tk::Dialog;

my $VERSION = 0.005;         # version of the Tk interface
my $version;                 # version of MELTING
my $blackhole;               # guess what ...
my $NNDIR;                   # address of the files containing the calorimetric data 

#----------------------
# Arguments of melting
#----------------------

my $altNN      = "";         # option -A  -> alternative set of calorimetric parameters
my $complement = "";         # option -C  -> complementary sequence
my $factor     = "4";        # option -F  -> factor correcting the nucleic acid concentration
my $magnesium  = "0.0";      # option -G  -> magnesium concentration
my $type       = "dnadna";   # option -H  -> type of hybridisation (dna/dna, dna/rna, rna/rna)
# option I useless here (infile)
my $Korrection = "san98a";   # option -K  -> salt correction
my $mismatches = "";         # option -M  -> alternative set of calorimetric parameters fot mismatches 
my $salt       = "0.0";      # option -N  -> salt concentration
# option O useless here? (outfile)
my $probe      = "0.0";      # option -P -> nuceic acid concentration
my $potassium  = "0.0";      # option -p -> potassium concentration
my $sequence   = "";         # option -S -> Sequence of one strand
my $tris       = "0.0";      # option -t -> tris concentration

#---------
# Results
#---------

my $enthalpy;              # enthalpy of the helix-coil transition
my $entropy;               # entropy of the helix-coil transition
my $tm;                    # melting-point (temperature of mid-fusion)

my $mw = MainWindow->new();
$mw->title("Melting");
$mw->resizable(0,0);
$version = @{ [split(" ",`melting -V`)] }[1];

###############################
# Presentation of the program #
###############################

my $presentFrame = $mw->Frame->pack(-expand => '1', 
				    -fill   => 'both'
				    );

$presentFrame->Label(-text       => "MELTING v$version
 Computation of the Tm according to the nearest-neighbour method 
        Copyright (C) Nicolas Le Nov�re 1997-2001",
		     -relief     => 'raised',
		     -background => '#DDDDDD',
		     )->pack(-expand => '1', 
			     -fill   => 'both'
			     );

######################################
# Enter the arguments of the program #
######################################

my $argFrame = $mw->Frame(-borderwidth => 1, 
			  -relief      => 'ridge'
			  )->pack(-expand => '1', 
				  -fill   => 'both'
				  );

#-----------
# sequences
#-----------

$argFrame->Label(-text => 'sequence => A,G,C,T,U,I,- <option -S>'
		 )->pack(-side   => 'top',
			 -expand => '1', 
			 -fill   => 'x'
			 );
$argFrame->Entry(-width        => '30',
		 -textvariable => \$sequence,
		 -font         => "-adobe-courier-medium-r-normal--0-0-0-0-m-0-iso8859-1",
		 -background   => '#999999',
		 -foreground   => '#EEEEEE'
		 )->pack(-side => 'top',
			 -expand => '1', 
			 -fill => 'x'
			 );

$argFrame->Label(-text => 'complement => A,G,C,T,U,I,- (only if mismatches, or inosine mismatches) <option -C>',
		 )->pack(-side   => 'top',
			 -expand => '1', 
			 -fill   => 'x'
			 );
$argFrame->Entry(-width        => '30',
		 -textvariable => \$complement,
		 -font         => "-adobe-courier-medium-r-normal--0-0-0-0-m-0-iso8859-1",
		 -background   => '#999999',
		 -foreground   => '#EEEEEE'
		 )->pack(-side   => 'top',
			 -expand => '1', 
			 -fill   => 'x'
			 );

my $numericalFrame = $argFrame->Frame->pack(-side   => 'right',
					    -expand => '1', 
					    -fill   => 'both'
					      );

#----------------------------
# nucleic acid concentration
#----------------------------

$numericalFrame->Label(-text => "[strand] <option -P>"
		   )->grid(-column => 0, 
			   -row    => 0,
			   ); 
$numericalFrame->Entry(-width        => 10,
		       -textvariable => \$probe,
		       -background   => '#999999',
		       -foreground   => '#EEEEEE'
		       )->grid(-column => 1, 
			       -row    => 0,
			       );
$numericalFrame->Label(-text => ']0..0.1] M.l-1'
		       )->grid(-column => 2, 
			       -row    => 0,
			       );

#--------------------
# ions concentration 
#--------------------

$numericalFrame->Label(-text => "[salt] <option -N>"
		  )->grid(-column => 0, 
			  -row    => 1,
			   );
$numericalFrame->Entry(-width        => 10,
		  -textvariable => \$salt,
		  -background   => '#999999',
		  -foreground   => '#EEEEEE'
		  )->grid(-column => 1, 
			  -row    => 1,
			  );
$numericalFrame->Label(-text => ']0...10] M.l-1'
		  )->grid(-column => 2, 
			  -row    => 1,
			  );
			  
$numericalFrame->Label(-text => "[potassium] <option -k>"
		  )->grid(-column => 0, 
			  -row    => 2,
			   );
$numericalFrame->Entry(-width        => 10,
		  -textvariable => \$potassium,
		  -background   => '#999999',
		  -foreground   => '#EEEEEE'
		  )->grid(-column => 1, 
			  -row    => 2,
			  );
$numericalFrame->Label(-text => ' M.l-1'
		  )->grid(-column => 2, 
			  -row    => 2,
			  );

$numericalFrame->Label(-text => "[tris] <option -t>"
		  )->grid(-column => 0, 
			  -row    => 3,
			   );
$numericalFrame->Entry(-width        => 10,
		  -textvariable => \$tris,
		  -background   => '#999999',
		  -foreground   => '#EEEEEE'
		  )->grid(-column => 1, 
			  -row    => 3,
			  );
$numericalFrame->Label(-text => ' M.l-1'
		  )->grid(-column => 2, 
			  -row    => 3,
			  );
$numericalFrame->Label(-text => "[magnesium] <option -G>"
		  )->grid(-column => 0, 
			  -row    => 4,
			   );
$numericalFrame->Entry(-width        => 10,
		  -textvariable => \$magnesium,
		  -background   => '#999999',
		  -foreground   => '#EEEEEE'
		  )->grid(-column => 1, 
			  -row    => 4,
			  );
$numericalFrame->Label(-text => ' M.l-1'
		  )->grid(-column => 2, 
			  -row    => 4,
			  );

#--------------------------------------------------------------
# stoichiometric correction for the nucleic acid concentration
#--------------------------------------------------------------

$numericalFrame->Label(-text => "Is one strand in excess (PCR ...)?"
		       )->grid(-column     => 0,
			       -columnspan => 2,
			       -row        => 5,
			       );

my $self_Lb = $numericalFrame->Label(-text => "Self-complementary sequence?"
				     )->grid(-column     => 0, 
					     -columnspan => 2,
					     -row        => 6,
					     );

my $self_Cb = $numericalFrame->Checkbutton(-onvalue  => 1,
					   -offvalue => 4,
					   -variable => \$factor,
					   )->grid(-column => 2, 
						   -row    => 5,
						   );

$numericalFrame->Checkbutton(-onvalue  => 2,
			     -offvalue => 4,
			     -variable => \$factor,
			     -command  => sub{ 
				 if ($self_Cb->cget(-indicatoron)){
				     $self_Cb->configure(-state       => 'disabled',
							 -indicatoron => 0);
				 } else {
				     $self_Cb->configure(-state       => 'normal',
							 -indicatoron => 1);
				 }
			     },
			     )->grid(-column => 2, 
				     -row    => 6,
				     );
#------------------
# Salt corrections
#------------------

my $saltCorrFrame = $numericalFrame->Frame->grid(-row        => 8,
						 -column     => 0,
						 -columnspan => 3
						   );

$saltCorrFrame->Label(-text => "[salt] correction <option -K> (if only Na+ are present)" 
		      )->pack(-side => 'left'
			      );

foreach (qw(san98a san96a wet91a owc08)){ # added explicit Owc salt correction (MP)
    $saltCorrFrame->Radiobutton(-text        => $_, 
				-value       => $_, 
				-variable    => \$Korrection,
				-selectcolor => '#999999'
				)->pack(
					);
}
#----------------------------------------------------
# list the available alternative calorimetric tables
#----------------------------------------------------

# Open the directory containing the calorimetric data
# remove the name of the files related to mismatches
# and in the future, controlled by the option -M and not -A

$NNDIR = `melting -p`;
chomp $NNDIR;
$NNDIR =~ s/^\s*path:\s+//;
$NNDIR =~ s/\/$//;

opendir NNF, $NNDIR or print "could not open the directory $NNDIR\n";

my @files = grep {($_ !~ /(mm|de)\.nn/) && ($_ !~ /(san05a|bre07a)\.nn/) && $_ !~ /^\./} readdir NNF;
close NNF;

my $AltSetMenu = $argFrame->Menubutton(-text   => "Alt NN sets\n<option -A>", 
				       -relief => 'raised'
				       )->pack(-side => 'bottom'
					       );

foreach my $file (@files){
    $AltSetMenu->command(-label   => $file,    
			 -command => sub{$altNN=$file;}
			 );
}
 
# Note that another possibility is the use of Tk::Optionmenu
# However, it does not exhibit any -text option (the visible value is the first
# of the list) of -label

#------------------
# hybridation type
#------------------

my $HybTypeFrame = $argFrame->Frame->pack(-expand => '1', 
					  -fill   => 'both',
					  -side   => 'left'
					    );

$HybTypeFrame->Label(-text => "Hybridation type\n<option -H>"
		     )->pack(
			     );

foreach (qw(DNA/DNA DNA/RNA RNA/RNA)){
    $HybTypeFrame->Radiobutton(-text        => $_, 
			       -value       => $_, 
			       -variable    => \$type,
			       -selectcolor => '#999999'
			       )->pack(
				       );
}

#####################################
# Output the results of the program #
#####################################

my $resFrame = $mw->Frame(-borderwidth => 1, 
			  -relief      => 'ridge'
			  )->pack(-side   => 'left',
				  -expand => '1', 
				  -fill   => 'both'
				  );

#--------------------------------
# The explanation of each result
#--------------------------------

$resFrame->Label(-text => 'enthalpy'
		    )->grid(-row    => 0,
			    -column => 0
			    );
$resFrame->Label(-text => 'entropy'
		    )->grid(-row    => 1,
			    -column => 0
			    );
$resFrame->Label(-text => 'melting temperature'
		    )->grid(-row    => 2,
			    -column => 0
			    );

#---------------------------------
# Finally, the results themselves
#---------------------------------

$resFrame->Entry(-width        => 10,
		 -textvariable => \$enthalpy
		 )->grid(-row    => 0,
			 -column => 1
			 );
$resFrame->Entry(-width        => 10,
		 -textvariable => \$entropy  
		 )->grid(-row    => 1,
			 -column => 1
			   );
$resFrame->Entry(-width        => 10,
		 -textvariable => \$tm
		 )->grid(-row    => 2,
			 -column => 1
			   );

#-------------------------
# The unit of each result
#-------------------------

$resFrame->Label(-text => 'J.mol-1'
		 )->grid(-row    => 0,
			 -column => 2
			 );
$resFrame->Label(-text => 'J.mol-1.K-1'
		 )->grid(-row    => 1,
			 -column => 2
			  );
$resFrame->Label(-text => '� C'
		 )->grid(-row    => 2,
			 -column => 2
			  );

##############
# Here we go #
##############

my $commandFrame = $mw->Frame(-borderwidth => 1, 
			      -relief      => 'ridge',
			      )->pack(-side   => 'left',
				      -expand => '1', 
				      -fill   => 'both'
				      );

$commandFrame->Button(-text    => 'Help',
		      -width   => '10',
		      -command => sub{show_help()}
		      )->pack(
			      );

$commandFrame->Button(-text    => 'run',
		      -width   => '10',
		      -command => sub{compute()}
		      )->pack(
			      );

$commandFrame->Button(-text    => "quit",
		      -width   => '10',
		      -command => sub { exit } 
		      )->pack(
			      );

# Note that tk exit function is used here instead of
# -command => [$main => 'destroy']
# The code after the MainLoop is therefore ignored.
# Remind it if you add code after the GUI termination

MainLoop;

#####################################################
# Open a toplevel window with the manpage of melting
#####################################################

sub show_help{
    if(!Exists(my $helpWindow)){
	$helpWindow = $mw->Toplevel();
	$helpWindow->title("Melting help");
	my $helpText = $helpWindow->Scrolled("Text", 
					     -font => "-adobe-courier-medium-r-normal--0-0-0-0-m-0-iso8859-1"
					     )->pack(-expand => 1,
						     -fill   => 'both'
						   );
	open(FH,"melting.hlp") || die "unable to open melting.hlp";
	while (<FH>){
	    $helpText->insert('end', $_);
	}
	close(FH);
	$helpWindow->Button(-text    => "close",
			    -command => sub{$helpWindow->withdraw}
			    )->pack(
				    );
    } else {
	$helpWindow->deiconify();
	$helpWindow-> raise();
    }
}

#########################################
# Run melting with the chosen parameters 
#########################################

sub compute {
    my $options;
    my @results;
    my $line;
    my $error;
    my @content;
    
#--------------------------------------
# Construction de la ligne de commande
#--------------------------------------
    $type =~ s!/!!;
    $type =~ tr/A-Z/a-z/;
    $options = "-H$type";
    if ($Korrection ne "owc08"){
    	$options.=" -K$Korrection";} # pass salt correction method to command line (MP)

    if ($altNN ne ""){
	$options.=" -A$altNN";
    }
    
    if ($sequence =~ /[^agctui-]/i or $sequence eq ""){ # allow inosine (MP)
	message("error","error","The sequence is empty or contains one or more illegal digits. The only legal digits are a,A,c,C,g,G,i,I,t,T,u,U");
	return;
    } else {$options.=" -S$sequence";}
    
    if ($complement ne ""){
	if ($complement =~ /[^agctui-]/i){ # allow inosine (MP)
	    message("error","error","The complementary sequence contains one or more illegal digits. The only legal digits are a,A,c,C,g,G,i,I,t,T,u,U");
	    return;
	} elsif (length $sequence != length $complement){
	    message("error","error","The sequence and its complement have different length");
	    return;
	} else {$options.=" -C$complement";}
    }

    if ($probe <= 0 or $probe > 0.1){
	message("error","error","The concentration of nucleic acid is out of range. It has to be over 0 and under 0.1 M");
	return;
    } else {$options.=" -P$probe";}
    
    if ($salt < 0 or $salt > 10){
	message("error","error","The concentration of salt out of range. It has to be over 0 and under 10 M");
	return;
    } else {$options.=" -N$salt";}
    
    if ($potassium < 0){
	message("error","error","The concentration of potassium must be positive.");
	return;
    } elsif ($potassium > 0) {$options.=" -k$potassium";} # Check for non-zero before passing to command line (MP)
    
    if ($tris < 0){
	message("error","error","The concentration of tris must be positive.");
	return;
    } elsif ($tris > 0){$options.=" -t$tris";}  # Check for non-zero before passing to command line (MP)
    
    if ($magnesium < 0){
	message("error","error","The concentration of magnesium must be positive.");
	return;
    } elsif ($magnesium > 0 or $Korrection eq "owc08") {$options.=" -G$magnesium";}  # Check for non-zero before passing to command line. Also force Owc by passing magnesium even if zero (MP)
    
    if ($factor ne "" and $factor ne "default"){
	if ($factor =~ /\D/){
	    message("error","error","The factor of correction for the nucleic acid concentration seems to contain one or more illegal characters. It has to be a number");
	    return;
	} else {$options.=" -F$factor";}
    }
    
#-------------
# C'est parti
#-------------
    
    @results=`melting $options -v -q 2>&1`;
    print "@results";

#------------------------
# Parsing of the results
#------------------------
    
    foreach $line (@results){
	if ($line =~ /^ \w/){
	    $error.=$line;
	} elsif ($line =~ /^  \w/){
	    @content = split(" ",$line);
	    if ($content[0] eq "Enthalpy:"){
		$enthalpy = $content[1];
	    } elsif ($content[0] eq "Entropy:"){
		$entropy = $content[1];
	    } elsif ($content[0] eq "Melting"){
		$tm = $content[2];	
	    } else{
		# Do nothing
	    }
	}
    }
    if (defined $error){
	$error =~ s/\n//;
	message("error","error",$error);
    }
}

#############################################################################
# This function opens a window and display the messages required by argument 
#############################################################################

sub message{

    my $title = shift;  # chops the title as the first argument
    my $bitmap = shift; # chops the bitmap as the second argument
    my $message = shift;
    
    my $message_window = $mw->Dialog(
				     -title => $title,
				     -bitmap => $bitmap,
				     -text  => $message,
				     -default_button => 'OK',
				     -buttons        => ['OK']
				     )->Show;  
}
