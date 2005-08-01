#!/usr/bin/perl -w

#***********************************************************************#
#
#     This program retrives a list of amplicon donor pairs with edited
#     phd files given a directory skeleton
# 
#***********************************************************************#

package ProcessingObjects::Metrics;

=pod

=head1 NAME

ProcessingObjects::Metrics

=head1 AUTHOR

Anushka Nethisinghe

=head1 DESCRIPTION

Object contains functions to calculate 
    percent amplicon coverage
    percent high minors in SNR clear range
    percent peak ratio average
    percent high minors overall
    number of bases in opposite strand

=cut

use Config::IniFiles;
use File::Path;
use File::Copy;
use File::Basename;
use File::Spec;
use IO::Handle;
use Getopt::Std;
use Date::Format;

use Data::Dumper;
###################################################################

sub new{
    my $class = shift;
    my $object = {};
    if(@_ == 5){

	# CHECK FOR ALL 5 REQUIRED INPUT PARAMETERS
	$object->{'gene'} = shift;
	$object->{'amplicon'} = shift;
	$object->{'donor'} = shift;
	$object->{'record'} = shift;
	$object->{'cfg_file'} = shift;
	
        # CHECK FOR VALID CONFIGURATION FILE
	if (!(-e $object->{'cfg_file'}  && -r $object->{'cfg_file'}) || (-d $object->{'cfg_file'})){
	    print STDERR "Config file ($opt_c) can not be read. Quitting . . . \n";
	    return $object;
	}
	
	# READ INI FILE
	my %config;
	tie %config, 'Config::IniFiles', ( -file => "$object->{cfg_file}" );

	# COPY NECESSARY CONFIG PARAMETERS TO METRICS OBJECT AND SET DIRECTORY VARIABLES
	$object->{'bl2seq'} = $config{'Blast_Tools'}{'bl2seq_exe'};
	
	my $target_study_dir = $config{'Target_Info'}{'target_dir_root'}."/".$config{'Target_Info'}{'target_project_type'}."/".$config{'Target_Info'}{'target_project'}."/".$config{'Target_Info'}{'target_study'};
	my $work_dir = $target_study_dir."/".$config{'Target_Info'}{'target_work'};
	my $data_dir = $target_study_dir."/".$config{'Target_Info'}{'target_data_dir_name'};
	my $primer_design_dir = $data_dir."/".$config{'PrimerDesign_Info'}{'primer_design_dir_name'};
	my $sequencing_dir = $data_dir."/".$config{'Target_Info'}{'target_sequencing_dir_name'};
	my $manifest_dir = $sequencing_dir."/".$config{'Target_Info'}{'target_manifest_dir_name'};
	my $manifest_files = $manifest_dir."/".$config{'Manifest_Info'}{'manifest_query_file'}.".txt";
	my @pda_list = split(",", $config{'PrimerDesign_Info'}{'source_dir_root'});

	my $pda_fasta = my $amplicons_dir = '';
	
	foreach my $dir (@pda_list){
	    $pda_fasta = `ls -1 $dir/*/*/$object->{'gene'}/$config{'PrimerDesign_Info'}{'reference_amplicons_dir_name'}/$object->{'amplicon'}.fasta`;
	    if($pda_fasta =~ /ls\s+No\s+match\./){
		next;
	    }else{
		chomp($pda_fasta);
		$pda_fasta =~ m/($dir)(\S+)($object->{'gene'})\//;
		$amplicons_dir = $dir.$2.$object->{'gene'};
		$reference_amplicons_dir = $amplicons_dir."/".$config{'PrimerDesign_Info'}{'reference_amplicons_dir_name'};
		$object->{'pda_no_primers_fasta'} = `ls -1 $dir/*/*/$object->{'gene'}/$config{'PrimerDesign_Info'}{'reference_amplicons_dir_name'}/$object->{'amplicon'}_no_primers.fasta | tail -1`;
		chomp($object->{'pda_no_primers_fasta'});
		$object->{'pda_no_primers_length'} = `$config{'Emboss_Tools'}{'infoseq_exe'} $object->{'pda_no_primers_fasta'} -only -length`;
		last;
	    }
	}
	my $reference_amplicons_dir = $amplicons_dir."/".$config{'PrimerDesign_Info'}{'reference_amplicons_dir_name'};
	my $reference_amplicons_blastdb = $reference_amplicons_dir ."/".$config{'PrimerDesign_Info'}{'all_amplicons_file'};
	my $primer_design_report_file = $amplicons_dir."/".$config{'PrimerDesign_Info'}{'primer_info_file'};
	my $primer_design_template_file = $amplicons_dir."/template.fasta";

	$curr_dir = $work_dir."/".$object->{'gene'}."/".$object->{'amplicon'}."/".$object->{'donor'};

	my @rec = split(/\t/, $object->{'record'});
	my $check_amp = $rec[3];
	my $check_dna = $rec[4];
	chomp($check_dna);
	chomp($check_amp);
	$check_amp =~ s/\s+//;
	$check_dna =~ s/\s+//;
	if($check_amp eq $object->{'amplicon'} && $check_dna eq $object->{'donor'}){
	    my $plate = $rec[0];
	    my $row = $rec[1];
	    my $col = $rec [2];

	    my $ztr_file = `ls -1 $sequencing_dir/*/*/$rec[0]_$rec[1]$rec[2]_*`;
	    chomp($ztr_file);
	    my $temp_name = basename($ztr_file);
	    my @ztr_split = split("_", $temp_name);
	    my $new_name = '';
	    if( $ztr_file =~ m/forward/){
		$object->{'ztr_file'} = $object->{'donor'}."_f_".$object->{'amplicon'}."_".$ztr_split[5].".ztr";
		$object->{'fc'} = 1;
		$object->{'psym'} = "TF";
	    }else{
		$object->{'ztr_file'} = $object->{'donor'}."_r_".$object->{'amplicon'}."_".$ztr_split[5].".ztr";
		$object->{'fc'} = 0;
		$object->{'psym'} = "TR";
	    }
	    my ($tracefile, $path, $type) = fileparse($object->{'ztr_file'}, ".ztr");
	    my $traceid = $ztr_split[5];
	    chomp($traceid);
            my $scfname = $tracefile.".scf";
	    my $snrfile = $config{'TraceAnalysisDir_Info'}{'snr'}."/".$scfname.".snr";
            my $phdfile = $config{'TraceAnalysisDir_Info'}{'phd'}."/".$scfname.".phd.1";
	    $object->{'polyfile'} = $config{'TraceAnalysisDir_Info'}{'poly'}."/".$scfname.".poly";
	    if (-e $snrfile){
		if (-e $phdfile){
		    if (-e $object->{'polyfile'}){
			$object->{'snr_s'} = `cut -d ' ' -f 3 $snrfile`;
			$object->{'snr_e'} = `cut -d ' ' -f 4 $snrfile`;
			$object->{'fastafile'} = $config{'TraceAnalysisDir_Info'}{'fasta'}."/".$scfname.".fasta.screen";
			$object->{'metricsfile'} = $config{'TraceAnalysisDir_Info'}{'metrics'}."/".$scfname.".metrics";
			$object->{'barcodeid'} = $plate;
			$object->{'traceid'} = $traceid;
			$object->{'wellrow'} = $row;
			$object->{'wellcol'} = $col;
			
		    }else{
			print STDERR "Polyfile, $polyfile, does not exist ($!)\n";
			return null;
		    }
		}else{
		    print STDERR "Phd file, $phdfile, does not exist ($!)\n";
		    return null;
		}
	    }else{
		print STDERR "SNR file, $snrfile, does not exist ($!)\n";
		return null;
	    }
	}else{
	    print STDERR "Invalid record $object->{'record'}) \n";
	    return null;
	}
    }else{
	print STDERR "Invalid inputs (@_)\n";
	return null;
    }
    bless $object, $class;
    $object->_initialize(@_);
    return $object;
}

###############################################################################

sub _initialize {
    my $object = shift;
}

###################################################################

#RETURN THE BARCODEID
sub getBarcodeId{
    my $metrics_object = shift;
    return $metrics_object->{'barcodeid'};
}

###################################################################

# RETURN THE TRACE ID
sub getTraceId{
    my $metrics_object = shift;
    return $metrics_object->{'traceid'};
}

###################################################################

# RETURN THE WELL COLUMN
sub getWellRow{
    my $metrics_object = shift;
    return $metrics_object->{'wellrow'};
}

###################################################################

# RETURN THE WELL ROW
sub getWellCol{
    my $metrics_object = shift;
    return $metrics_object->{'wellcol'};
}

###################################################################

# RETURN THE AMPLICON NAME
sub getAmplicon{
    my $metrics_object = shift;
    return $metrics_object->{'amplicon'};
}

###################################################################

# RETURN THE DONOR NAME
sub getDonor{
    my $metrics_object = shift;
    return $metrics_object->{'donor'};
}

###################################################################

# RETURN THE PRIMER SYMBOL (DIRECTION)
sub getPSYMDirection{
    my $metrics_object = shift;
    return $metrics_object->{'psym'};
}

###################################################################

# CALCULATE THE PERCENT OF THE AMPLICON THAT IS COVERED    
sub getPercentAmpliconCoverage{
    my $metrics_object = shift; 

    my $percentAmpliconCoverage = 0;
    my $num_hsps = `$metrics_object->{'bl2seq'} -p blastn -F F -e 1e-10 -D 1 -i $metrics_object->{'pda_no_primers_fasta'} -j $metrics_object->{'fastafile'} | grep -v "^#" | wc -l`;
    my @bl2seq_records = `$metrics_object->{'bl2seq'} -p blastn -F F -e 1e-10 -D 1 -i $metrics_object->{'pda_no_primers_fasta'} -j $metrics_object->{'fastafile'} | grep -v "^#"`;

    if ($num_hsps > 0){
	my $bl2seq = $direction = "";
	my $s = $e = $fs = $fe = $percent = 0;
	if ($metrics_object->{'fc'} > 0 ){ 
	    foreach $bl2seq_record (@bl2seq_records){
		@bl2seq_rec_split = split("\t", $bl2seq_record);
	   
		if((($bl2seq_rec_split[7]-$bl2seq_rec_split[6])*($bl2seq_rec_split[9]-$bl2seq_rec_split[8])) > 1){
		    if($s == 0){
			$s = $bl2seq_rec_split[8];
			$e = $bl2seq_rec_split[9];
		    }else{
			if($s>$bl2seq_rec_split[8]){
			    $s = $bl2seq_rec_split[8];
			}
			if($e<$bl2seq_rec_split[9]){
			    $e = $bl2seq_rec_split[9];
			}
		    }
		}
	    }
	    if(($s < $metrics_object->{'snr_e'}) && ($e > $metrics_object->{'snr_s'})){
		$fs = $s;
		if($fs < $metrics_object->{'snr_s'}){
		    $fs = $metrics_object->{'snr_s'};
		}
		$fe = $e;
		if($fe > $metrics_object->{'snr_e'}){
		    $fe = $metrics_object->{'snr_e'};
		}
	    }
	    $percent = (($fe-$fs)+1)*100/$metrics_object->{'pda_no_primers_length'};
	}else{
	    foreach $bl2seq_record (@bl2seq_records){
		@bl2seq_rec_split = split("\t", $bl2seq_record);

		if((($bl2seq_rec_split[7]-$bl2seq_rec_split[6])*($bl2seq_rec_split[9]-$bl2seq_rec_split[8])) < 1){
		    if($s == 0){
			$s = $bl2seq_rec_split[9];
			$e = $bl2seq_rec_split[8];
		    }else{
			if($s>$bl2seq_rec_split[9]){
			    $s = $bl2seq_rec_split[9];
			}
			if($e<$bl2seq_rec_split[8]){
			    $e = $bl2seq_rec_split[8];
			}
		    }
		}
	    }
	    if(($s < $metrics_object->{'snr_e'}) && ($e > $metrics_object->{'snr_s'})){
		$fs = $s;
		if($fs < $metrics_object->{'snr_s'}){
		    $fs = $metrics_object->{'snr_s'};
		}
		$fe = $e;
		if($fe > $metrics_object->{'snr_e'}){
		    $fe = $metrics_object->{'snr_e'};
		}
	    }
	    $percent = (($fe-$fs)+1)*100/$metrics_object->{'pda_no_primers_length'};
	}
	$percentAmpliconCoverage = sprintf("%.2f", $percent);
	return $percentAmpliconCoverage;
    }else{
	return 0;
    }
    
}
###################################################################

sub getPercentHighMinorsOverall{
    my $metrics_object = shift;
    
    my $percentHighMinorsOverall = 0;
    open(POLYHD, "$metrics_object->{'polyfile'}") || die "Cannot open $metrics_object->{'polyfile'} ($!)\n";
    my $maja = 1.0;
    my $rec_num = $a = 0;
    foreach $line (<POLYHD>){
	chomp($line);
	$rec_num++;
	my @line_split = split(" ", $line);
	$maja = 1.0;
	if($maja < $line_split[2]){
	    $maja = $line_split[2];
	}
	if($rec_num > 1){
	    if($line_split[0] ne $line_split[4]){
		if($line_split[4] ne "N"){
		    if(($line_split[6]/$maja) >0.2){
			$a++;
		    }
		}
	    }
	}
    }
    if($a > $metrics_object->{'pda_no_primers_length'}){
	$a = $metrics_object->{'pda_no_primers_length'};
    }
    $percentHighMinorsOverall = sprintf("%.2f", 100 * $a / $metrics_object->{'pda_no_primers_length'});
    close POLYHD;
    return $percentHighMinorsOverall;
}

###################################################################

# CALCULATE THE PERCENT OF HIGH MINORS IN THE SNR CLEAR RANGE
sub getPercentHighMinorsInSNRCLR{
    my $metrics_object = shift;
    
    my $percentHighMinorsInSNRCLR = 0;
    my $maja = 1.0;
    my $rec_num = $a = 0;
    open(POLYHD, "$metrics_object->{'polyfile'}") || die "Cannot open $metrics_object->{'polyfile'} ($!)\n";
    foreach $line (<POLYHD>){
	chomp($line);
	$rec_num++;
	my @line_split = split(" ", $line);
	$maja = 1.0;
	if($maja < $line_split[2]){
	    $maja = $line_split[2];
	}
	if($rec_num > 1){
	    if(($rec_num - 1) >= $metrics_object->{'snr_s'}){
		if(($rec_num - 1) <= $metrics_object->{'snr_e'}){
		    if($line_split[0] ne $line_split[4]){
			if($line_split[4] ne "N"){
			    if(($line_split[6]/$maja) >0.2){
				$a++;
			    }
			}
		    }
		}
	    }
	}
    }
    if($metrics_object->{'snr_s'} == $metrics_object->{'snr_e'} || $a > $metrics_object->{'pda_no_primers_length'}){
	$a = $metrics_object->{'pda_no_primers_length'};
    }
    $percentHighMinorsInSNRCLR = sprintf("%.2f", 100 * $a / $metrics_object->{'pda_no_primers_length'});
    close POLYHD;
    return $percentHighMinorsInSNRCLR;
}

###################################################################

# CALCULATE THE AVERAGE PEAK RATIO FOR A TRACE FILE
sub getPercentPeakRatioAverage{
    my $metrics_object = shift;
    
    my $percentPeakRatioAverage = 0;
    open(POLYHD, "$metrics_object->{'polyfile'}") || die "Cannot open $metrics_object->{'polyfile'} ($!)\n";
    my $maja = 1.0;
    my $rec_num = $nm = $rs = 0;
    foreach $line (<POLYHD>){
	chomp($line);
	$rec_num++;
	my @line_split = split(" ", $line);
	$maja = 1.0;
	if($maja < $line_split[2]){
	    $maja = $line_split[2];
	}
	if($rec_num > 1){
	    if(($rec_num - 1) >= $metrics_object->{'snr_s'}){
		if(($rec_num - 1) <= $metrics_object->{'snr_e'}){
		    if($line_split[0] ne $line_split[4]){
			if($line_split[4] ne "N"){
			    $nm++;
			    $rs += $line_split[6]/$maja
			}
		    }
		}
	    }
	}
    }
    if($nm == 0){
	$nm = 1;
    }
    $percentPeakRatioAverage = sprintf("%.2f", 100 * $rs / $nm);
    close POLYHD;
    return $percentPeakRatioAverage;
}

###################################################################

# CALCULATE THE NUMBER OF BASES OF OPPOSITE STRAND PRIMING
sub getNumBasesOppositeStrand{
    my $metrics_object = shift;
    
    my $numBasesOppositeStrand = 0;
    open(POLYHD, "$metrics_object->{'polyfile'}") || die "Cannot open $metrics_object->{'polyfile'} ($!)\n";
    $majFile = $metrics_object->{'polyfile'}.".fasta_major";
    $minFile = $metrics_object->{'polyfile'}.".fasta_minor";
    open(MAJOROUTHD, ">$majFile") || die "Cannot open $majFile ($!)\n";
    open(MINOROUTHD, ">$minFile") || die "Cannot open $minFile ($!)\n";
    my $rec_num = $fe = $fs = $s = $e = $nob = 0;

    # CREATE A FASTA FILE FOR THE MAJOR AND MINOR BASES
    foreach $line (<POLYHD>){
	chomp($line);
	$rec_num++;
	my @line_split = split(" ", $line);
	if($rec_num == 1){
	    printf MAJOROUTHD ">major_%s\n", $line_split[0];
	    printf MINOROUTHD ">minor_%s\n", $line_split[0];
	}else{
	    printf MAJOROUTHD "%s", $line_split[0];
	    if($line_split[4] eq "N"){
		printf MINOROUTHD "%s", $line_split[0];
	    }else{
		printf MINOROUTHD "%s", $line_split[4];
	    }
	}
    }

    close POLYHD;
    close MAJOROUTHD;
    close MINOROUTHD;

    my @bl2seq_records = `$metrics_object->{'bl2seq'} -p blastn -e 1e-10 -D 1 -i $majFile -j $minFile | grep -v "^#"`;

    foreach $bl2seq_record (@bl2seq_records){
	my @bl2seq_rec_split = split("\t", $bl2seq_record);

	# COUNT NUMBER OF OPPOSITE STRAND PRIMING BASES
	if((($bl2seq_rec_split[7]-$bl2seq_rec_split[6])*($bl2seq_rec_split[9]-$bl2seq_rec_split[8])) < 1){
	    if($s == 0){
		$s = $bl2seq_rec_split[9];
		$e = $bl2seq_rec_split[8];
	    }else{
		if($s>$bl2seq_rec_split[9]){
		    $s = $bl2seq_rec_split[9];
		}
		if($e<$bl2seq_rec_split[8]){
		    $e = $bl2seq_rec_split[8];
		}
	    }
	}
    }
    $fs = $s;
    $fe = $e;
    $nob = (($fe-$fs) + 1);
    if($s == 0){
	$nob = 0;
    }
    $numBasesOppositeStrand = $nob;
    return $numBasesOppositeStrand;
}


1;
