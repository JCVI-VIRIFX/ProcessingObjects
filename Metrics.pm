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

use strict;
use Config::IniFiles;
use File::Path;
use File::Basename;
use MLDBM 'DB_File';
use File::Spec;
use IO::Handle;
use Getopt::Std;
use Date::Format;

###################################################################

sub new{
    my $class = shift;
    my $object = {};
    if(@_ == 4){

	# CHECK FOR ALL 5 REQUIRED INPUT PARAMETERS
        $object->{'ztr_file'} = shift;
	$object->{'record'} = shift;
	$object->{'cfg_file'} = shift;
	$object->{'curr_work_dir'} = shift;

	chdir($object->{'curr_work_dir'});
	
	# Set core variables
	my @rec = split(/\t/, $object->{'record'});
	$object->{'barcodeid'} = $rec[0];
	$object->{'wellrow'} = $rec[1];
	$object->{'wellcol'} = $rec[2];
	$object->{'amplicon'} = $rec[3];
	$object->{'dna'} = $rec[4];
	$object->{'psym'} = $rec[5];
	
	my @ztr_split = split("_", $object->{'ztr_file'});
	if($object->{'ztr_file'} =~ /\_(\d+).ztr/){
	    $object->{'traceid'} = $1;
	}
	if($object->{'psym'} eq "TF"){
	    $object->{'fc'} = 1;
	}else{
	    $object->{'fc'} = 0;
	}
	my ($tracefile, $path, $type) = fileparse($object->{'ztr_file'}, ".ztr");

        # Check for valid configuration file
	if (!(-e $object->{'cfg_file'}  && -r $object->{'cfg_file'}) || (-d $object->{'cfg_file'})){
	    print STDERR "Config file ($object->{'cfg_file'}) can not be read. Quitting . . . \n";
	    return $object;
	}
	
	# Read ini file
	my %config;
	tie %config, 'Config::IniFiles', ( -file => "$object->{cfg_file}" );

	# Set directory and exe variables
	$object->{'bl2seq'} = $config{'Blast_Tools'}{'bl2seq_exe'};
	$object->{'blast'} = $config{'Blast_Tools'}{'blast_exe'};
	$object->{'repeat_detector'} = $config{'AmpliconHandling_Tools'}{'repeat_detector_exe'};
	$object->{'fuzznuc'} = $config{'Emboss_Tools'}{'fuzznuc_exe'};
	$object->{'fasta_dir'} = $config{'TraceAnalysisDir_Info'}{'fasta'};
	$object->{'all_amplicons_file'} = $config{'PrimerDesign_Info'}{'all_amplicons_file'};
	$object->{'blast_dir'} = $config{'TraceAnalysisDir_Info'}{'blast'};
	
	# Retrieve primer design reference amplicon data path form list of primer design data paths
	my $source_dir_root = $config{'PrimerDesign_Info'}{'source_dir_root'};
	my %all_amp_hash = ();
	tie(%all_amp_hash, 'MLDBM',$source_dir_root) || die "Cannot to to $source_dir_root ($!)\n";
	my $pda_fasta = $all_amp_hash{$object->{'amplicon'}};
	my ($temp_fasta_file, $reference_amplicons_dir) = fileparse($pda_fasta, ".fasta");
	my $amplicons_dir = $reference_amplicons_dir;
	$amplicons_dir =~ s/reference_amplicons\///g;
	$object->{'pda_no_primers_fasta'} = $reference_amplicons_dir."/".$object->{'amplicon'}."_no_primers.fasta";
	$object->{'pda_no_primers_length'} = `$config{'Emboss_Tools'}{'infoseq_exe'} $object->{'pda_no_primers_fasta'} -only -length`;
	$object->{'ref_amp_blastdb'} = $reference_amplicons_dir ."/".$config{'PrimerDesign_Info'}{'all_amplicons_file'};
	$object->{'primer_info_file'} = $amplicons_dir."/".$config{'PrimerDesign_Info'}{'primer_info_file'};
	my $primer_design_template_file = $amplicons_dir."/template.fasta";

	
	my $scfname = $tracefile.".scf";
	my $snrfile = $config{'TraceAnalysisDir_Info'}{'snr'}."/".$scfname.".snr";
	my $phdfile = $config{'TraceAnalysisDir_Info'}{'phd'}."/".$scfname.".phd.1";
	$object->{'polyfile'} = $config{'TraceAnalysisDir_Info'}{'poly'}."/".$scfname.".poly";
	$object->{'fuzznuc_fasta'} = "$object->{'curr_work_dir'}/$config{'TraceAnalysisDir_Info'}{'fasta'}/$scfname.fasta";
	if (-e $snrfile){
	    if (-e $phdfile){
	        if (-e $object->{'polyfile'}){
		    $object->{'snr_s'} = `cut -d ' ' -f 3 $snrfile`;
		    $object->{'snr_e'} = `cut -d ' ' -f 4 $snrfile`;
		    $object->{'fastafile'} = $config{'TraceAnalysisDir_Info'}{'fasta'}."/".$scfname.".fasta.screen";
		    $object->{'metricsfile'} = $config{'TraceAnalysisDir_Info'}{'metrics'}."/".$scfname.".metrics";
		}else{
		    print STDERR "Polyfile, $object->{'polyfile'}, does not exist ($!)\n";
		    return 1;
		}
	    }else{
		print STDERR "Phd file, $phdfile, does not exist ($!)\n";
		return 1;
	    }
	}else{
	    print STDERR "SNR file, $snrfile, does not exist ($!)\n";
	    return 1;
	}
    }else{
	print STDERR "Invalid inputs (@_)\n";
	return 1;
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

# CALCULATE AND RETURN THE MATCHING REFERENCE SEQUENCE
sub getMatchingAmpliconReferenceSequence{
    my $metrics_object = shift;

    chdir($metrics_object->{'fasta_dir'});

    my $scf_screen;
    my $scf_screen = `ls *$metrics_object->{'traceid'}.scf.fasta.screen`;
    chomp($scf_screen);
    my($file, $path, $type) = fileparse($scf_screen, '.scf.fasta.screen');
    chomp($file);

    my $blast_out = "../blast_dir/".$file.".blastn.out";
    my @hit_count = `grep "^>$metrics_object->{'amplicon'}" $blast_out`;

    # IF AMPLICON MATCHES REFERENCE CREATE SCREEN AND QUAL FILES
    if ($#hit_count >= 0){
	my $strand = `grep "Strand" $blast_out | head -1`;
	if($strand =~ m/Strand\s=\s\S+\s\/\s(\S+)/){
	    my $direction = $1;
	    if($direction eq "Plus" &&  $metrics_object->{'psym'} eq "TR"){
		$metrics_object->{'hit_amp'} = "reference(TF)";
	    }elsif($direction eq "Minus" && $metrics_object->{'psym'} eq "TF"){
		$metrics_object->{'hit_amp'} = "reference(TR)";
	    }else{
		$metrics_object->{'hit_amp'} = "-8";
	    }
	}
    }else{

	@hit_count = `grep "^>" $blast_out`;
	
	# IF AMPLICON MATCHES ANOTHER REFERENCE SEQUENCE IN THE SAME DIRECTORY THEN TAG AMPLICON
	if($#hit_count >= 0){
	    
	    # GET STRAND DIRECTION
	    my @strand = `grep "Strand" $blast_out`;
	    $strand[0] =~ m/Strand\s=\s\S+\s\/\s(\S+)/;
	    my $direction = $1;
	    
	    # GET AMPLICON HIT NAME
	    $hit_count[0] =~ />(\S+)_reference\.scf/;
	    my $amp_match = $1;
	    
	    if($direction eq "Plus"){
		$metrics_object->{'hit_amp'} = $amp_match."(TF)";
	    }else{
		$metrics_object->{'hit_amp'} = $amp_match."(TR)";
	    }
	}else{
	    
	    # ELSE CHECK IF AMPLICON MATCHES REFERENCE OF AMPLICONS IN THE SAME PROJECT
	    $metrics_object->{'ref_amp_blastdb'} =~ m/(\S+)(high|strict)/;
	    my $proj_dir = $1;
	    my $ref_proj_blastdb = $proj_dir.$metrics_object->{'all_amplicons_file'};
	    my $blast_all_project = "../".$metrics_object->{'blast_dir'}."/".$file.".blastn_all_project.out";

	    `$metrics_object->{'blast'} -p blastn -d $ref_proj_blastdb -i $scf_screen -o $blast_all_project -e 1e-10 -F F -v 5 -b 5`;
	    @hit_count = `grep "^>" $blast_all_project`;
	    
	    
	    if($#hit_count >= 0){
		
		# GET STRAND DIRECTION
		my @strand = `grep "Strand" $blast_all_project`;
		$strand[0] =~ m/Strand\s=\s\S+\s\/\s(\S+)/;
		my $direction = $1;
		
		# GET AMPLICON HIT NAME
		$hit_count[0] =~ />(\S+)_reference\.scf/;
		my $amp_match = $1;
		
		if($direction eq "Plus"){
		    $metrics_object->{'hit_amp'} = $amp_match."(TF)";
		}else{
		    $metrics_object->{'hit_amp'} = $amp_match."(TR)";
		}
	    }else{
		my $ref_all_blastdb = $proj_dir."../../".$metrics_object->{'all_amplicons_file'};
		my $blast_all = "../".$metrics_object->{'blast_dir'}."/".$file.".blastn_all.out";

		`$metrics_object->{'blast'} -p blastn -d $ref_all_blastdb -i $scf_screen -o $blast_all -e 1e-10 -F F -v 5 -b 5`;
		@hit_count = `grep "^>" $blast_all`;
		if($#hit_count >= 0){
		    
		    # GET STRAND DIRECTION

		    my @strand = `grep "Strand" $blast_all`;
		    $strand[0] =~ m/Strand\s=\s\S+\s\/\s(\S+)/;
		    my $direction = $1;
		    
		    # GET AMPLICON HIT NAME
		    $hit_count[0] =~ />(\S+)_reference\.scf/;
		    my $amp_match = $1;
		    
		    if($direction eq "Plus"){
			$metrics_object->{'hit_amp'} = $amp_match."(TF)";
		    }else{
			$metrics_object->{'hit_amp'} = $amp_match."(TR)";
		    }
		}else{
		    $metrics_object->{'hit_amp'} = "NO_MATCH";
		}
	    }
	}
    }
    chdir(File::Spec->updir);
    return $metrics_object->{'hit_amp'};
}

###################################################################

# CALCULATE THE PERCENT OF THE AMPLICON THAT IS COVERED    
sub getPercentAmpliconCoverage{
    my $metrics_object = shift; 

    my @ref_start_stutter = ();
    my $ref_start_stutter = 0;
    my $i;

   # DETERMINE START STUTTER IN REFERENCE SEQUENCE
    @ref_start_stutter = `$metrics_object->{'repeat_detector'} -f $metrics_object->{'pda_no_primers_fasta'} -s 13 -m 22 | grep "^[1-9]" | cut -f 1`;
    for($i=0; $i<=$#ref_start_stutter; $i++){
	my $ref_start_stut = $ref_start_stutter[$i];
	chomp($ref_start_stut);
	if($ref_start_stut >10){
	    $ref_start_stutter = $ref_start_stut;
	    $i = $#ref_start_stutter+1 ;
	}
    }

    my $percentAmpliconCoverage = 0;
    my $bases_for_full_coverage = $metrics_object->{'pda_no_primers_length'};
    if ( $ref_start_stutter > 0 ) {
      $bases_for_full_coverage = $ref_start_stutter;
    }

    my $num_hsps = `$metrics_object->{'bl2seq'} -p blastn -F F -e 1e-10 -D 1 -i $metrics_object->{'pda_no_primers_fasta'} -j $metrics_object->{'fastafile'} | grep -v "^\#" | wc -l`;
    my @bl2seq_records = `$metrics_object->{'bl2seq'} -p blastn -F F -e 1e-10 -D 1 -i $metrics_object->{'pda_no_primers_fasta'} -j $metrics_object->{'fastafile'} | grep -v "^\#"`;

    if ($num_hsps > 0){
	my($bl2seq, $direction) = "";
	my($s, $e, $fs, $fe, $percent) = 0;
	if ($metrics_object->{'fc'} > 0 ){ 
	    foreach my $bl2seq_record (@bl2seq_records){
		my @bl2seq_rec_split = split("\t", $bl2seq_record);
	   
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
	    $percent = (($fe-$fs)+1)*100/$bases_for_full_coverage;
	}else{
	    foreach my $bl2seq_record (@bl2seq_records){
		my @bl2seq_rec_split = split("\t", $bl2seq_record);

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
	    $percent = (($fe-$fs)+1)*100/$bases_for_full_coverage;
	}
        if ( $percent > 100 ) {
            $percent = 100.0;
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
    foreach my $line (<POLYHD>){
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
    foreach my $line (<POLYHD>){
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
    my($rec_num, $nm, $rs) = 0;
    foreach my $line (<POLYHD>){
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
    my $fuzznuc = "/usr/local/packages/EMBOSS-2.9.0/bin/fuzznuc";
    my $fuzznuc_cmd = "";
    my $numBasesOppositeStrand = my $max_base = my $hit_max_base = 0;
    my $min_amp_length = 100;
    my $max_amp_length = 850;
    my $seq_primer_pigtail_length = 19;

    # Retrieve primer sequences
    my $fwd_primer = `grep $metrics_object->{'amplicon'} $metrics_object->{'primer_info_file'} | cut -f 6`;
    chomp($fwd_primer);
    $fwd_primer =~ tr/[AaCcGgTtMmRrYyKkVvHhDdBb]/[TtGgCcAaKkYyRrMmBbDdHhVv]/;
    my $rev_comp_fwd_primer = reverse($fwd_primer);
    my $rev_primer = `grep $metrics_object->{'amplicon'} $metrics_object->{'primer_info_file'} | cut -f 7`;
    chomp($rev_primer);
    $rev_primer =~ tr/[AaCcGgTtMmRrYyKkVvHhDdBb]/[TtGgCcAaKkYyRrMmBbDdHhVv]/;
    my $rev_comp_rev_primer = reverse($rev_primer);

    # Determine end of valid sequence in trace
    if($metrics_object->{'psym'} eq "TF"){
	$fuzznuc_cmd = "$metrics_object->{'fuzznuc'}
				 -sequence $metrics_object->{'fuzznuc_fasta'}
                                 -pattern $rev_comp_rev_primer
                                 -mismatch 2
                                 -outfile $metrics_object->{'fuzznuc_fasta'}.fuzznuc";
    }else{
	$fuzznuc_cmd = "$metrics_object->{'fuzznuc'} 
				 -sequence $metrics_object->{'fuzznuc_fasta'}
                                 -pattern $rev_comp_fwd_primer
                                 -mismatch 2
                                 -outfile $metrics_object->{'fuzznuc_fasta'}.fuzznuc";
    }
    $fuzznuc_cmd =~ s/\n/ /g;
    `$fuzznuc_cmd`;

    my @base_info_arr = `grep -v "^#" $metrics_object->{'fuzznuc_fasta'}.fuzznuc | grep -v "Start"`;
    foreach my $base_info (@base_info_arr){
	$base_info =~ s/(\s+)/:/g;
	my @max_base_info_arr = split(":", $base_info);
	my $seq_length = $max_base_info_arr[1];
	if($seq_length > $min_amp_length){
	    $max_base = $seq_length + $seq_primer_pigtail_length;
	    $hit_max_base = 1;
	}
    }

    if(!($hit_max_base)){
	$max_base = $max_amp_length;
    }

    open(POLYHD, "$metrics_object->{'polyfile'}") || die "Cannot open $metrics_object->{'polyfile'} ($!)\n";
    my $majFile = $metrics_object->{'polyfile'}.".fasta_major";
    my $minFile = $metrics_object->{'polyfile'}.".fasta_minor";
    open(MAJOROUTHD, ">$majFile") || die "Cannot open $majFile ($!)\n";
    open(MINOROUTHD, ">$minFile") || die "Cannot open $minFile ($!)\n";
    my($rec_num, $fe, $fs, $s, $e, $nob) = 0;

    # Generate major and minor sequence fasta files
    foreach my $line (<POLYHD>){
	chomp($line);
	$rec_num++;
	my @line_split = split(" ", $line);
	if($rec_num == 1){
	    printf MAJOROUTHD ">major_%s\n", $line_split[0];
	    printf MINOROUTHD ">minor_%s\n", $line_split[0];
	}elsif($rec_num <= $max_base+1){
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

    foreach my $bl2seq_record (@bl2seq_records){
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
