#!/usr/bin/perl -w

#***********************************************************************#
#
#     For a given study XML file, parses the data for accessing 
# 
#***********************************************************************#

package ProcessingObjects::ReadStudyXML;

=pod

=head1 NAME

ProcessingObjects::ReadStudyXML;

=head1 AUTHOR

Anushka Brownley

=head1 DESCRIPTION

This module handles the parsing of a study XML file and creates and object
to access the data.
	
The following is a list of the subroutines in the module:
=cut

use strict;
use XML::Simple;
use Getopt::Std;

use Data::Dumper;
###################################################################

sub new{
    my $class = shift;
    my $object = {};
    if(@_ == 1){

	# Store input parameters
	$object->{'xml_file'} = shift;
    }else{
	die "Invalid inputs (@_)\n";
    }
    bless $object, $class;
    if($object->_initialize(@_)){
	die "ERROR Initializing object\n";
    }
    
    return $object;
}

###############################################################################

# Initialize object
sub _initialize {
    my $object = shift;

    # Initialize analysis variables
    $object->{'xml_hash'};
    $object->{'curr_template_count'} = 0;
    $object->{'curr_amplicon_count'} = 0;
    $object->{'curr_dna_count'} = 0;
    $object->{'total_dna_count'} = 0;

    if($object->_loadXML()){
	print STDERR "ERROR loading XML file\n";
	return 1;
    }
    return 0;
}

###################################################################
# Loading XML data
sub _loadXML{
		
=head2 int _loadXML();

	This function load the XML data into a ReadStudyXML object
=cut
    my $object = shift;

    if(!(-e $object->{'xml_file'})){
	print STDERR "ERROR XML file does not exist\n";
	return 1;
    }
    
    my $xmlObj = new XML::Simple;
    $object->{'xml_hash'} = $xmlObj->XMLin($object->{'xml_file'}, KeyAttr=>'name', ContentKey => '-content');
    
    return 0;
}

###################################################################
# Returns the study name 
sub getStudy{
		
=head2 string getStudy();

	This function returns the study name
=cut
    my $object = shift;

    return $object->{'xml_hash'}->{'name'};
}

###################################################################
# Return 1 or 0 if amplicon is or is not a valid amplicon for 
# this study
sub isValidAmplicon{
		
=head2 int isValidAmplicon(string ampliconName);

	This function returns 1 or 0 if amplicon is or is not a 
	valid amplicon for this study
=cut
    my $object = shift;
    my $amp_name = shift;
    
    if(defined $object->{'xml_hash'}->{'ampliconList'}->{'amplicon'}->{$amp_name} ||
       $object->{'xml_hash'}->{'ampliconList'}->{'amplicon'}->{'name'} eq $amp_name){
	return 1;
    }else{
	return 0;
    }
}

###################################################################
# Return list of all amplicons
sub getAllAmplicons{
		
=head2 arrRef getAllAmplicons();

	This function returns a reference to an array of all amplicons
=cut
    my $object = shift;

    my @amp_arr;
    if(defined $object->{'xml_hash'}->{'ampliconList'}->{'amplicon'}->{'name'}){
	push @amp_arr, $object->{'xml_hash'}->{'ampliconList'}->{'amplicon'}->{'name'};
    }else{
        @amp_arr = keys(%{$object->{'xml_hash'}->{'ampliconList'}->{'amplicon'}});
	if($#amp_arr < 1){
	    print STDERR "ERROR: returning amplicon list\n";
	    return 1;
	}
    }
    
    return \@amp_arr;
}

###################################################################
# Return list of stutter amplicons
sub getStutterAmplicons{
		
=head2 arrRef getStutterAmplicons(in StutterLengthCutoff);

	This function returns a reference to an array of stutter amplicons
	with a stutter length greater than or equal to the input cutoff
=cut
    my $object = shift;
    my $stutter_cutoff = shift;
    
    if($stutter_cutoff <0 || $stutter_cutoff ne int($stutter_cutoff)){
	print STDOUT "Please enter a valid integer for stutter cutoff\n";
	return 1;
    }
    
    print "Retrieving stutter amplicons with cutoff $stutter_cutoff\n";
    my @amp_arr;
    if(defined $object->{'xml_hash'}->{'ampliconList'}->{'amplicon'}->{'name'}){
	if($object->{'xml_hash'}->{'ampliconList'}->{'amplicon'}->{'stutter'} >= $stutter_cutoff){
	    push @amp_arr, $object->{'xml_hash'}->{'ampliconList'}->{'amplicon'}->{'name'};
	}
    }else{
	foreach my $amp (keys(%{$object->{'xml_hash'}->{'ampliconList'}->{'amplicon'}})){
	    if($object->{'xml_hash'}->{'ampliconList'}->{'amplicon'}->{$amp}->{'stutter'} >= $stutter_cutoff){
		    push @amp_arr, $amp;
	    }
	}
    }
    
    return \@amp_arr;
}

###################################################################
# Return the status of a given amplicon
sub getAmpliconStatus{
		
=head2 string getAmpliconStatus(string AmpliconName);

	This function returns the status of a given amplicon
=cut
    my $object = shift;
    my $amp_name = shift;

    if(defined $object->{'xml_hash'}->{'ampliconList'}->{'amplicon'}->{'name'}){
	if($object->{'xml_hash'}->{'ampliconList'}->{'amplicon'}->{'name'} eq $amp_name){
	    return $object->{'xml_hash'}->{'ampliconList'}->{'amplicon'}->{'status'};
	}else{
	    print "ERROR: $amp_name does not exist as amplicon name\n";
	    return 1;
	}
    }else{
	if(defined($object->{'xml_hash'}->{'ampliconList'}->{'amplicon'}->{$amp_name})){
	    return $object->{'xml_hash'}->{'ampliconList'}->{'amplicon'}->{$amp_name}->{'status'};
	}else{
	    print "ERROR: $amp_name does not exist as amplicon name\n";
	    return 1;
	}
    }
    
    return 1;
}

###################################################################
# Return the true or false whether given amplicon is a stutter
# amplicon based on the cutoff
sub getAmpliconStutterStatus{
		
=head2 string getAmpliconStutterStatus(string AmpliconName, int stutterCutoff);

	This function returns true or false whether given amplicon 
	is a stutter amplicon based on the input cutoff
=cut
    my $object = shift;
    my $amp_name = shift;
    my $stutter_cutoff = shift;

    if($stutter_cutoff <0 || $stutter_cutoff ne int($stutter_cutoff)){
	print STDOUT "Please enter a valid integer for stutter cutoff\n";
	return 1;
    }
    
    if(defined $object->{'xml_hash'}->{'ampliconList'}->{'amplicon'}->{'name'}){
	if($object->{'xml_hash'}->{'ampliconList'}->{'amplicon'}->{'name'} eq $amp_name){
	    if($object->{'xml_hash'}->{'ampliconList'}->{'amplicon'}->{'stutter'} >= $stutter_cutoff){
		return 'true';
	    }else{
		return 'false';
	    }
	}else{
	    print "ERROR: $amp_name does not exist as amplicon name\n";
	    return 1;
	}
    }else{
	if(defined($object->{'xml_hash'}->{'ampliconList'}->{'amplicon'}->{$amp_name})){
	    if($object->{'xml_hash'}->{'ampliconList'}->{'amplicon'}->{$amp_name}->{'stutter'} >= $stutter_cutoff){
		return 'true';
	    }else{
		return 'false';
	    }
	}else{
	    print "ERROR: $amp_name does not exist as amplicon name\n";
	    return 1;
	}
    }
    
    return 1;
}

###################################################################
# Return 1 or 0 if dna is or is not a valid dna for 
# this study
sub isValidDNA{
		
=head2 int isValidDNA(string dnaName);

	This function returns 1 or 0 if dna is or is not a 
	valid dna for this study
=cut
    my $object = shift;
    my $dna_name = shift;
    
    if(defined $object->{'xml_hash'}->{'dnaList'}->{'dna'}->{$dna_name} ||
       $object->{'xml_hash'}->{'dnaList'}->{'dna'}->{'name'} eq $dna_name){
	return 1;
    }else{
	return 0;
    }
}

###################################################################
# Return list of all dnas
sub getAllDNAs{
		
=head2 arrRef getAllDNAs();

	This function returns a reference to an array of all DNAs
=cut
    my $object = shift;

    my @dna_arr;
    if(defined $object->{'xml_hash'}->{'dnaList'}->{'dna'}->{'name'}){
	push @dna_arr, $object->{'xml_hash'}->{'dnaList'}->{'dna'}->{'name'};
    }else{
        @dna_arr = keys(%{$object->{'xml_hash'}->{'dnaList'}->{'dna'}});
	if($#dna_arr < 1){
	    print STDERR "ERROR: returning DNA list\n";
	    return 1;
	}
    }
    
    return \@dna_arr;
}

###################################################################
# Return the status of a given dna
sub getDNAStatus{
		
=head2 string getDNAStatus(string DNAName);

	This function returns the status of a given dna
=cut
    my $object = shift;
    my $dna_name = shift;

    my @dna_arr;
    if(defined $object->{'xml_hash'}->{'dnaList'}->{'dna'}->{'name'}){
	if($object->{'xml_hash'}->{'dnaList'}->{'dna'}->{'name'} eq $dna_name){
	    return $object->{'xml_hash'}->{'dnaList'}->{'dna'}->{'status'};
	}else{
	    print "ERROR: $dna_name does not exist as dna name\n";
	    return 1;
	}
    }else{
	if(defined($object->{'xml_hash'}->{'dnaList'}->{'dna'}->{$dna_name})){
	    return $object->{'xml_hash'}->{'dnaList'}->{'dna'}->{$dna_name}->{'status'};
	}else{
	    print "ERROR: $dna_name does not exist as dna name\n";
	    return 1;
	}
    }
    
    return 1;
}

###################################################################
# Returns true or false if the input array does or does not match
# the complete set of DNAs for the study
sub checkDNASetComplete{
		
=head2 string checkDNASetComplete(arrRef dnaArr);

	Returns true or false if the input array does or does 
	not match the complete set of DNAs for the study
=cut
    my $object = shift;
    my $in_dna_arr_ref = shift;
    
    my @in_dna_arr = @{$in_dna_arr_ref};
    if($#in_dna_arr < 0){
	print STDERR "ERROR: input DNA array is empty\n";
	return 1;
    }
    
    my $dna_arr_ref = $object->getAllDNAs();
    my @dna_arr = @{$dna_arr_ref};
    my @sort_dna_arr = sort {$a <=> $b} @dna_arr;
    my @sort_in_dna_arr = sort {$a <=> $b} @in_dna_arr;
    if($#sort_dna_arr != $#sort_in_dna_arr){
	return 'false';
    }
    for(my $i=0; $i<=$#sort_dna_arr; $i++){
	if($sort_dna_arr[$i] ne $sort_in_dna_arr[$i]){
	    return 'false';
	}
    }
    
    return 'true';
}

###################################################################
# Return 1 or 0 if template is or is not a valid template for 
# this study
sub isValidTemplate{
		
=head2 int isValidTemplate(string templateName);

	This function returns 1 or 0 if template is or is not a 
	valid template for this study
=cut
    my $object = shift;
    my $template_name = shift;
    
    if(defined $object->{'xml_hash'}->{'templateList'}->{'template'}->{$template_name} ||
       $object->{'xml_hash'}->{'templateList'}->{'template'}->{'name'} eq $template_name){
	return 1;
    }else{
	return 0;
    }
}

###################################################################
# Return list of all templates
sub getAllTemplates{
		
=head2 arrRef getAllTemplates();

	This function returns a reference to an array of all templates
=cut
    my $object = shift;

    my @template_arr;
    if(defined $object->{'xml_hash'}->{'templateList'}->{'template'}->{'name'}){
	push @template_arr, $object->{'xml_hash'}->{'templateList'}->{'template'}->{'name'};
    }else{
        @template_arr = keys(%{$object->{'xml_hash'}->{'templateList'}->{'template'}});
	if ($#template_arr < 1){
	    print STDERR "ERROR: returning template list\n";
	    return 1;
	}
    }
    
    return \@template_arr;
}

###################################################################
# Return true or false for conducting variant manual review
sub conductVariantReview{
		
=head2 string conductVariantReview();

	This function returns true or false for conducting variant manual review
=cut
    my $object = shift;

    if(defined $object->{'xml_hash'}->{'manualReview'}->{'variantReview'}){
        return $object->{'xml_hash'}->{'manualReview'}->{'variantReview'}->{'status'};
    }else{
	print STDERR "ERROR: variant review in XML file does not exist\n";
	return 1;
    }
}

###################################################################
# Return true or false for conducting indel manual review
sub conductIndelReview{
		
=head2 string conductIndelReview();

	This function returns true or false for conducting indel manual review
=cut
    my $object = shift;

    if(defined $object->{'xml_hash'}->{'manualReview'}->{'indelReview'}){
        return $object->{'xml_hash'}->{'manualReview'}->{'indelReview'}->{'status'};
    }else{
	print STDERR "ERROR: indel review in XML file does not exist\n";
	return 1;
    }
}

###################################################################
# Return true or false for conducting stutter manual review
sub conductStutterReview{
		
=head2 string conductStutterReview();

	This function returns true or false for conducting stutter manual review
=cut
    my $object = shift;

    if(defined $object->{'xml_hash'}->{'manualReview'}->{'stutterReview'}){
        return $object->{'xml_hash'}->{'manualReview'}->{'stutterReview'}->{'status'};
    }else{
	print STDERR "ERROR: stutter review in XML file does not exist\n";
	return 1;
    }
}

###################################################################
# Return true or false for conducting somatic manual review
sub conductSomaticReview{
		
=head2 string conductSomaticReview();

	This function returns true or false for conducting somatic manual review
=cut
    my $object = shift;

    if(defined $object->{'xml_hash'}->{'manualReview'}->{'somaticReview'}){
        return $object->{'xml_hash'}->{'manualReview'}->{'somaticReview'}->{'status'};
    }else{
	print STDERR "ERROR: somatic review in XML file does not exist\n";
	return 1;
    }
}

###################################################################
# Return true or false for conducting delivery manual review
sub conductDataDeliveryReview{
		
=head2 string conductDataDeliveryReview();

	This function returns true or false for conducting delivery manual review
=cut
    my $object = shift;

    if(defined $object->{'xml_hash'}->{'manualReview'}->{'deliveryReview'}){
        return $object->{'xml_hash'}->{'manualReview'}->{'deliveryReview'}->{'status'};
    }else{
	print STDERR "ERROR: delivery review in XML file does not exist\n";
	return 1;
    }
}

###################################################################
# Return nav file for conducting variant manual review
sub variantReviewNavFile{
		
=head2 string variantReviewNavFile();

	This function returns nav file for conducting variant manual review
=cut
    my $object = shift;

    if($object->{'xml_hash'}->{'manualReview'}->{'variantReview'}->{'status'} eq 'true'){
        return $object->{'xml_hash'}->{'manualReview'}->{'variantReview'}->{'reviewExtent'};
    }else{
	print STDERR "Variant review should not take place for this study\n";
	return 1;
    }
}

###################################################################
# Return nav file for conducting indel manual review
sub indelReviewNavFile{
		
=head2 string indelReviewNavFile();

	This function returns nav file for conducting indel manual review
=cut
    my $object = shift;

    if($object->{'xml_hash'}->{'manualReview'}->{'indelReview'}->{'status'} eq 'true'){
        return $object->{'xml_hash'}->{'manualReview'}->{'indelReview'}->{'reviewExtent'};
    }else{
	print STDERR "Indel review should not take place for this study\n";
	return 1;
    }
}

###################################################################
# Return nav file for conducting stutter manual review
sub stutterReviewNavFile{
		
=head2 string stutterReviewNavFile();

	This function returns nav file for conducting stutter manual review
=cut
    my $object = shift;

    if($object->{'xml_hash'}->{'manualReview'}->{'stutterReview'}->{'status'} eq 'true'){
        return $object->{'xml_hash'}->{'manualReview'}->{'stutterReview'}->{'reviewExtent'};
    }else{
	print STDERR "Stutter review should not take place for this study\n";
	return 1;
    }
}

###################################################################
# This is the destructor function
sub DESTROY{
		
=head2  DESTROY();

	This is the destructor function
=cut
}

1;
