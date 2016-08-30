package PatternModule::TFBSAnno;

use strict;
use warnings;
use List::Util qw(max);

=head2 createVariantTFBSHash

  Args       : None
  Example    : my %variantInformation = createVariantTFBSHash()
  Description: Creates a hash that stores information on the variant's effects on 
               transcription factor binding sites
  Returntype : Hash

=cut

sub createVariantTFBSHash {
	my %data = ('ENSTFBSANNO_TFBS1_NAME'						=> "NA",
				'ENSTFBSANNO_TFBS1_MATRIX'						=> "NA",
				'ENSTFBSANNO_TFBS1_POS'							=> ".",
				'ENSTFBSANNO_TFBS1_LENGTH' 						=> ".",
				'ENSTFBSANNO_TFBS1_DIST_TO_START' 				=> ".",
				'ENSTFBSANNO_TFBS1_SIZE_OF_INDEL' 				=> ".",
				'ENSTFBSANNO_TFBS1_REF_SCORE' 					=> ".",
				'ENSTFBSANNO_TFBS1_ALT_SCORE' 					=> ".",
				'ENSTFBSANNO_TFBS1_INFO_CONTENT'				=> "NA",
				'ENSTFBSANNO_TFBS1_CONS_SCORE' 					=> ".",
				'ENSTFBSANNO_TFBS2_NAME'						=> "NA",
				'ENSTFBSANNO_TFBS2_MATRIX'						=> "NA",
				'ENSTFBSANNO_TFBS2_POS'							=> ".",
				'ENSTFBSANNO_TFBS2_LENGTH' 						=> ".",
				'ENSTFBSANNO_TFBS2_DIST_TO_START' 				=> ".",
				'ENSTFBSANNO_TFBS2_SIZE_OF_INDEL' 				=> ".",
				'ENSTFBSANNO_TFBS2_REF_SCORE' 					=> ".",
				'ENSTFBSANNO_TFBS2_ALT_SCORE' 					=> ".",
				'ENSTFBSANNO_TFBS2_INFO_CONTENT'				=> "NA",
				'ENSTFBSANNO_TFBS2_CONS_SCORE' 					=> ".",
				'ENSTFBSANNO_TFBS3_NAME'						=> "NA",
				'ENSTFBSANNO_TFBS3_MATRIX'						=> "NA",
				'ENSTFBSANNO_TFBS3_POS'							=> ".",
				'ENSTFBSANNO_TFBS3_LENGTH' 						=> ".",
				'ENSTFBSANNO_TFBS3_DIST_TO_START' 				=> ".",
				'ENSTFBSANNO_TFBS3_SIZE_OF_INDEL' 				=> ".",
				'ENSTFBSANNO_TFBS3_REF_SCORE' 					=> ".",
				'ENSTFBSANNO_TFBS3_ALT_SCORE' 					=> ".",
				'ENSTFBSANNO_TFBS3_INFO_CONTENT'				=> "NA",
				'ENSTFBSANNO_TFBS3_CONS_SCORE' 					=> ".",
				'ENSTFBSANNO_TFBS4_NAME'						=> "NA",
				'ENSTFBSANNO_TFBS4_MATRIX'						=> "NA",
				'ENSTFBSANNO_TFBS4_POS'							=> ".",
				'ENSTFBSANNO_TFBS4_LENGTH' 						=> ".",
				'ENSTFBSANNO_TFBS4_DIST_TO_START' 				=> ".",
				'ENSTFBSANNO_TFBS4_SIZE_OF_INDEL' 				=> ".",
				'ENSTFBSANNO_TFBS4_REF_SCORE' 					=> ".",
				'ENSTFBSANNO_TFBS4_ALT_SCORE' 					=> ".",
				'ENSTFBSANNO_TFBS4_INFO_CONTENT'				=> "NA",
				'ENSTFBSANNO_TFBS4_CONS_SCORE' 					=> ".",
				'ENSTFBSANNO_TFBS5_NAME'						=> "NA",
				'ENSTFBSANNO_TFBS5_MATRIX'						=> "NA",
				'ENSTFBSANNO_TFBS5_POS'							=> ".",
				'ENSTFBSANNO_TFBS5_LENGTH' 						=> ".",
				'ENSTFBSANNO_TFBS5_DIST_TO_START' 				=> ".",
				'ENSTFBSANNO_TFBS5_SIZE_OF_INDEL' 				=> ".",
				'ENSTFBSANNO_TFBS5_REF_SCORE' 					=> ".",
				'ENSTFBSANNO_TFBS5_ALT_SCORE' 					=> ".",
				'ENSTFBSANNO_TFBS5_INFO_CONTENT'				=> "NA",
				'ENSTFBSANNO_TFBS5_CONS_SCORE' 					=> ".",
				'ENSTFBSANNO_TFBS6_NAME'						=> "NA",
				'ENSTFBSANNO_TFBS6_MATRIX'						=> "NA",
				'ENSTFBSANNO_TFBS6_POS'							=> ".",
				'ENSTFBSANNO_TFBS6_LENGTH' 						=> ".",
				'ENSTFBSANNO_TFBS6_DIST_TO_START' 				=> ".",
				'ENSTFBSANNO_TFBS6_SIZE_OF_INDEL' 				=> ".",
				'ENSTFBSANNO_TFBS6_REF_SCORE' 					=> ".",
				'ENSTFBSANNO_TFBS6_ALT_SCORE' 					=> ".",
				'ENSTFBSANNO_TFBS6_INFO_CONTENT'				=> "NA",
				'ENSTFBSANNO_TFBS6_CONS_SCORE' 					=> ".",
				'ENSTFBSANNO_TFBS7_NAME'						=> "NA",
				'ENSTFBSANNO_TFBS7_MATRIX'						=> "NA",
				'ENSTFBSANNO_TFBS7_POS'							=> ".",
				'ENSTFBSANNO_TFBS7_LENGTH' 						=> ".",
				'ENSTFBSANNO_TFBS7_DIST_TO_START' 				=> ".",
				'ENSTFBSANNO_TFBS7_SIZE_OF_INDEL' 				=> ".",
				'ENSTFBSANNO_TFBS7_REF_SCORE' 					=> ".",
				'ENSTFBSANNO_TFBS7_ALT_SCORE' 					=> ".",
				'ENSTFBSANNO_TFBS7_INFO_CONTENT'				=> "NA",
				'ENSTFBSANNO_TFBS7_CONS_SCORE' 					=> ".",
				'ENSTFBSANNO_TFBS8_NAME'						=> "NA",
				'ENSTFBSANNO_TFBS8_MATRIX'						=> "NA",
				'ENSTFBSANNO_TFBS8_POS'							=> ".",
				'ENSTFBSANNO_TFBS8_LENGTH' 						=> ".",
				'ENSTFBSANNO_TFBS8_DIST_TO_START' 				=> ".",
				'ENSTFBSANNO_TFBS8_SIZE_OF_INDEL' 				=> ".",
				'ENSTFBSANNO_TFBS8_REF_SCORE' 					=> ".",
				'ENSTFBSANNO_TFBS8_ALT_SCORE' 					=> ".",
				'ENSTFBSANNO_TFBS8_INFO_CONTENT'				=> "NA",
				'ENSTFBSANNO_TFBS8_CONS_SCORE' 					=> ".",
				'ENSTFBSANNO_TFBS9_NAME'						=> "NA",
				'ENSTFBSANNO_TFBS9_MATRIX'						=> "NA",
				'ENSTFBSANNO_TFBS9_POS'							=> ".",
				'ENSTFBSANNO_TFBS9_LENGTH' 						=> ".",
				'ENSTFBSANNO_TFBS9_DIST_TO_START' 				=> ".",
				'ENSTFBSANNO_TFBS9_SIZE_OF_INDEL' 				=> ".",
				'ENSTFBSANNO_TFBS9_REF_SCORE' 					=> ".",
				'ENSTFBSANNO_TFBS9_ALT_SCORE' 					=> ".",
				'ENSTFBSANNO_TFBS9_INFO_CONTENT'				=> "NA",
				'ENSTFBSANNO_TFBS9_CONS_SCORE' 					=> ".",
				'ENSTFBSANNO_TFBS10_NAME'						=> "NA",
				'ENSTFBSANNO_TFBS10_MATRIX'						=> "NA",
				'ENSTFBSANNO_TFBS10_POS'						=> ".",
				'ENSTFBSANNO_TFBS10_LENGTH' 					=> ".",
				'ENSTFBSANNO_TFBS10_DIST_TO_START' 				=> ".",
				'ENSTFBSANNO_TFBS10_SIZE_OF_INDEL' 				=> ".",
				'ENSTFBSANNO_TFBS10_REF_SCORE' 					=> ".",
				'ENSTFBSANNO_TFBS10_ALT_SCORE' 					=> ".",
				'ENSTFBSANNO_TFBS10_INFO_CONTENT'				=> "NA",
				'ENSTFBSANNO_TFBS10_CONS_SCORE' 				=> ".",
	);	
	return %data;
}

=head2 fetchVariantTFBSinfo

  Args       : "VariationFeature", "MotifFeature Adaptor", "Slice Adaptor"
  Example    : my $string = fetchVariantTFBSinfo($vf, $motifFeature_adaptor, $slice_adaptor)
  Description: Creates a string with information on what transcription factor binding sites 
  			   overlap a variant, the conservation score of the TFBS and how the predicted 
  			   binding is affected by the variant. The predicted binding is calculated both
  			   for the reference sequence and the alternative sequence. The scoring of the
  			   reference sequence is straight forward, as it is just scored using the 
  			   relative_affinity method. For alternative sequences all adjoining sequences 
  			   are, however, scanned for an alternative binding site. The site with the best
  			   predicted binding is then chosen as the transcription factor binding site of
  			   the alternative allele. For SNVs a total of motif-length + 2 sequences are 
  			   scanned, for insertions a total of motif-length + insertion-length + 1
  			   sequences are scanned and for deletions a total of motif-length + 1 sequences
  			   are scanned.
  Returntype : String

=cut

sub fetchVariantTFBSinfo {
	my ($vf, $motifFeature_adaptor, $slice_adaptor, $csscore_adaptor, $mlss) = @_;
	my @allelicString = split("/", $vf->allele_string());
	my $variantSlice = $slice_adaptor->fetch_by_region('chromosome',$vf->slice->seq_region_name(),$vf->start(),$vf->start()-1+length($allelicString[0]));
	my @motif_features = @{$motifFeature_adaptor->fetch_all_by_Slice($variantSlice)};
	my @keys = ();
	
	my $annotationString = ();
	
	# Getting information on transcription factor binding modules that are hit by variants
	if (scalar @motif_features > 0) {
		my %motifScoreHash = ();
		my %motifScoreSortingHash = ();
	
			for my $motif (@motif_features) {
				if ($vf->start() <= $motif->seq_region_end() && $motif->seq_region_start() <= $vf->end()) {
					my $slice = $slice_adaptor->fetch_by_region('chromosome', $vf->slice->seq_region_name(), $motif->seq_region_start(), $motif->seq_region_end());				
									
					# For SNPs
					if ($vf->allele_string() !~ /-/ && length($vf->allele_string()) == 3) {
						if ($motif->strand() == 1) {
							# Calculating binding scores
							my $refString = $slice->seq();
							my $refScore = $motif->binding_matrix->relative_affinity($refString);
							
							# Creating a sequence including the SNV
							my $SNVSlice = $slice_adaptor->fetch_by_region('chromosome', $vf->slice->seq_region_name(), $motif->seq_region_start() - 1000, $motif->seq_region_end() + 1000);
							my $SNVString = $SNVSlice->seq();
							substr($SNVString, $vf->start() - $motif->seq_region_start() + 1000, 1) = substr($vf->allele_string(), 2, 1);

							my @newBindingScores = ();	
											
							# Scoring all adjoining sequences and taking the highest scoring sequence as the new TFBS 					
							for (my $i = 0; $i < $motif->length + 2; $i++) {
								my $altString = substr($SNVString, 1000 + $vf->start() - $motif->seq_region_start() - $motif->length() + $i, $motif->length());
								push @newBindingScores, $motif->binding_matrix->relative_affinity($altString);
							}
							my $totalScore = max(@newBindingScores)-$refScore;
							
							# Motif start
							my $StartPos = $motif->seq_region_start();
							
							# Length of motif
							my $motiflength = $motif->binding_matrix->length();
						
							# Distance from variant to start of motif
							my $DistToStart = $vf->start() - $motif->seq_region_start();
						
							# Size of indel
							my $SizeOfIndel = "NA";
					
							# Name of motif
							my @motifName = $motif->display_label();

							# Conservation score calculation
							my $scores = $csscore_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $motiflength);
							my ($consScore) = calculateConservationScore($scores);
							
							# TFBS information content
							my $infoCont = ();
							for (my $i = 0; $i < $motif->length(); $i++) {
								if ($i != $motif->length() - 1) {
									$infoCont .= $motif->binding_matrix->{'ic'}->[$i] . ",";
								}
								if ($i == $motif->length() - 1) {
									$infoCont .= $motif->binding_matrix->{'ic'}->[$i];
								}
							}
						
							# Placing information into Hash. If the same motif has two scores for a variant, the most severe score is stored
							unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]}{Difference} <= $totalScore) {
								$motifScoreHash{$motifName[0]}{Start} = $StartPos;
								$motifScoreHash{$motifName[0]}{Length} = $motiflength;
								$motifScoreHash{$motifName[0]}{DistToStart} = $DistToStart;
								$motifScoreHash{$motifName[0]}{SizeOfIndel} = $SizeOfIndel;
								$motifScoreHash{$motifName[0]}{Original} = $refScore;
								$motifScoreHash{$motifName[0]}{New} = max(@newBindingScores);
								$motifScoreHash{$motifName[0]}{Difference} = $totalScore;
								$motifScoreHash{$motifName[0]}{InfoContent} = $infoCont;
								$motifScoreHash{$motifName[0]}{ConsScore} = $consScore;
							}
						}
				
						if ($motif->strand() == -1) {
							# Calculating binding scores
							my $refString = PatternModule::RegAnno::createComplementarySequence($slice->seq());
							my $refScore = $motif->binding_matrix->relative_affinity($refString);
							
							# Creating a sequence including the SNV
							my $SNVSlice = $slice_adaptor->fetch_by_region('chromosome', $vf->slice->seq_region_name(), $motif->seq_region_start() - 1000, $motif->seq_region_end() + 1000);
							my $SNVString = $SNVSlice->seq();
							substr($SNVString, $vf->start() - $motif->seq_region_start() + 1000, 1) = substr($vf->allele_string(), 2, 1);
							my @newBindingScores = ();	
											
							# Scoring all adjoining sequences and taking the highest scoring sequence as the new TFBS 					
							for (my $i = 0; $i < $motif->length + 2; $i++) {
								my $altString = PatternModule::RegAnno::createComplementarySequence(substr($SNVString, 1000 + $vf->start() - $motif->seq_region_start() - $motif->length() + $i, $motif->length()));
								push @newBindingScores, $motif->binding_matrix->relative_affinity($altString);
							}
							my $totalScore = max(@newBindingScores)-$refScore;	
							
							# Motif start
							my $StartPos = $motif->seq_region_start();
							
							# Length of motif
							my $motiflength = $motif->binding_matrix->length();
					
							# Distance from variant to start of motif
							my $DistToStart = $motif->seq_region_end() - $vf->start();
						
							# Size of indel
							my $SizeOfIndel = "NA";
						
							# Name of motif
							my @motifName = $motif->display_label();
							
							# Conservation score calculation
							my $scores = $csscore_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $motiflength);
							my ($consScore) = calculateConservationScore($scores);
							
							# TFBS information content
							my $infoCont = ();
							for (my $i = 0; $i < $motif->length(); $i++) {
								if ($i != $motif->length() - 1) {
									$infoCont .= $motif->binding_matrix->{'ic'}->[$i] . ",";
								}
								if ($i == $motif->length() - 1) {
									$infoCont .= $motif->binding_matrix->{'ic'}->[$i];
								}
							}
						
							unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]}{Difference} <= $totalScore) {
								$motifScoreHash{$motifName[0]}{Start} = $StartPos;
								$motifScoreHash{$motifName[0]}{Length} = $motiflength;
								$motifScoreHash{$motifName[0]}{DistToStart} = $DistToStart;
								$motifScoreHash{$motifName[0]}{SizeOfIndel} = $SizeOfIndel;
								$motifScoreHash{$motifName[0]}{Original} = $refScore;
								$motifScoreHash{$motifName[0]}{New} = max(@newBindingScores);
								$motifScoreHash{$motifName[0]}{Difference} = $totalScore;
								$motifScoreHash{$motifName[0]}{InfoContent} = $infoCont;
								$motifScoreHash{$motifName[0]}{ConsScore} = $consScore;
							}
						}	
					}
	
					# For Substitutions
					if ($vf->allele_string() !~ /-/ && length($vf->allele_string()) > 3) {
						if ($motif->strand() == 1) {
					
							# Calculating binding scores
							my $refString = $slice->seq();
							my $refScore = $motif->binding_matrix->relative_affinity($refString);
			
							# Creating a sequence including the substitution	
							# A slice is made with 1000 bases in each direction of the motif
							my $deletionSlice = $slice_adaptor->fetch_by_region('chromosome', $vf->slice->seq_region_name(), $motif->seq_region_start() - 1000, $motif->seq_region_end() + 1000);
							my $deletionString = $deletionSlice->seq();
							my @varrefalt = split("/", $vf->allele_string());
							substr($deletionString, $vf->start() - $motif->seq_region_start() + 1000, $vf->end() - $vf->start() + 1) = $varrefalt[1];
							my @newBindingScores = ();						
							for (my $i = 0; $i < $motif->length + 1; $i++) {
								my $altString = substr($deletionString, 1000 + $vf->start() - $motif->seq_region_start() - $motif->length() + $i, $motif->length());
								push @newBindingScores, $motif->binding_matrix->relative_affinity($altString);
							}
							my $totalScore = max(@newBindingScores)-$refScore;
						
							# Motif start
							my $StartPos = $motif->seq_region_start();
						
							# Length of motif
							my $motiflength = $motif->binding_matrix->length();
		
							# Distance from variant to start of motif
							my $DistToStart = -1;
							if ($vf->start() - $motif->seq_region_start() > -1) {
								$DistToStart = $vf->start() - $motif->seq_region_start();
							}
					
							# Size of substitution
							my @splitallele = split("/", $vf->allele_string());
							my $SizeOfsubstitution = length($splitallele[0]);
					
							# If the substitution is of the entire transcription factor binding motif
							if ($vf->start() - $motif->seq_region_start() <= 0 && $motif->seq_region_end() - $vf->end() <= 0) {
								$SizeOfsubstitution = $motif->binding_matrix->length();
							}
					
							# If the substitution start/end both lie within the transcription factor binding motif
							if ($vf->start() - $motif->seq_region_start() >= 0 && $motif->seq_region_end() - $vf->end() >= 0) {
								$SizeOfsubstitution = length($splitallele[0]);
							}
					
							# If the substitution start lies before the motif start, and the substitution end lies within the motif
							if ($vf->start() - $motif->seq_region_start() <= 0 && $motif->seq_region_end() - $vf->end() >= 0) {
								$SizeOfsubstitution = $vf->end() - $motif->seq_region_start() + 1;
							}
					
							# If the substitution start lies within motif, and the substitution end lies outside motif
							if ($vf->start() - $motif->seq_region_start() >= 0 && $motif->seq_region_end() - $vf->end() <= 0) {
								$SizeOfsubstitution = $motif->seq_region_end() - $vf->start() + 1;
							}								
					
							# Name of motif
							my @motifName = $motif->display_label();
						
							# Conservation score calculation
							my $scores = $csscore_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $motiflength);
							my ($consScore) = calculateConservationScore($scores);
						
							# TFBS information content
							my $infoCont = ();
							for (my $i = 0; $i < $motif->length(); $i++) {
								if ($i != $motif->length() - 1) {
									$infoCont .= $motif->binding_matrix->{'ic'}->[$i] . ",";
								}
								if ($i == $motif->length() - 1) {
									$infoCont .= $motif->binding_matrix->{'ic'}->[$i];
								}
							}
					
							unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]}{Difference} <= $totalScore) {
								$motifScoreHash{$motifName[0]}{Start} = $StartPos;
								$motifScoreHash{$motifName[0]}{Length} = $motiflength;
								$motifScoreHash{$motifName[0]}{DistToStart} = $DistToStart;
								$motifScoreHash{$motifName[0]}{SizeOfsubstitution} = $SizeOfsubstitution;
								$motifScoreHash{$motifName[0]}{Original} = $refScore;
								$motifScoreHash{$motifName[0]}{New} = max(@newBindingScores);
								$motifScoreHash{$motifName[0]}{Difference} = $totalScore;
								$motifScoreHash{$motifName[0]}{InfoContent} = $infoCont;
								$motifScoreHash{$motifName[0]}{ConsScore} = $consScore;
							}
						}
						if ($motif->strand == -1) {
					
							# Calculating binding scores
							my $refString = PatternModule::RegAnno::createComplementarySequence($slice->seq());
							my $refScore = $motif->binding_matrix->relative_affinity($refString);
			
							# Creating a sequence including the substitution
							# A slice is made with 1000 bases in each direction of the motif
							my $deletionSlice = $slice_adaptor->fetch_by_region('chromosome', $vf->slice->seq_region_name(), $motif->seq_region_start() - 1000, $motif->seq_region_end() + 1000);
							my $deletionString = $deletionSlice->seq();
							my @varrefalt = split("/", $vf->allele_string());
							substr($deletionString, $vf->start - $motif->seq_region_start() + 1000, $vf->end() - $vf->start() + 1) = $varrefalt[1];
							my @newBindingScores = ();						
							for (my $i = 0; $i < $motif->length() + 1; $i++) {
								my $altString = PatternModule::RegAnno::createComplementarySequence(substr($deletionString, 1000 + $vf->start() - $motif->seq_region_start() - $motif->length() + $i, $motif->length()));
								push @newBindingScores, $motif->binding_matrix->relative_affinity($altString);
							}
							my $totalScore = max(@newBindingScores)-$refScore;
						
							# Motif start
							my $StartPos = $motif->seq_region_start();
						
							# Length of motif
							my $motiflength = $motif->binding_matrix->length();

							# Distance from variant to start of motif
							my $DistToStart = -1;
							if ($motif->seq_region_end() - $vf->end() > -1) {
								$DistToStart = $motif->seq_region_end() - $vf->end();
							}
					
							# Size of substitution
							my @splitallele = split("/", $vf->allele_string());
							my $SizeOfsubstitution = length($splitallele[0]);
					
							# If the substitution is of the entire transcription factor binding motif
							if ($vf->start() - $motif->seq_region_start() <= 0 && $motif->seq_region_end() - $vf->end() <= 0) {
								$SizeOfsubstitution = $motif->binding_matrix->length();
							}
					
							# If the substitution start/end both lie within the transcription factor binding motif
							if ($vf->start() - $motif->seq_region_start() >= 0 && $motif->seq_region_end() - $vf->end() >= 0) {
								$SizeOfsubstitution = length($splitallele[0]);
							}
					
							# If the substitution start lies before the motif start, and the substitution end lies within the motif
							if ($vf->start() - $motif->seq_region_start() <= 0 && $motif->seq_region_end() - $vf->end() >= 0) {
								$SizeOfsubstitution = $vf->end() - $motif->seq_region_start() + 1;
							}
					
							# If the substitution start lies within motif, and the substitution end lies outside motif
							if ($vf->start() - $motif->seq_region_start() >= 0 && $motif->seq_region_end() - $vf->end() <= 0) {
								$SizeOfsubstitution = $motif->seq_region_end() - $vf->start() + 1;
							}		
					
							# Name of motif
							my @motifName = $motif->display_label();
						
							# Conservation score calculation
							my $scores = $csscore_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $motiflength);
							my ($consScore) = calculateConservationScore($scores);
						
							# TFBS information content
							my $infoCont = ();
							for (my $i = 0; $i < $motif->length(); $i++) {
								if ($i != $motif->length() - 1) {
									$infoCont .= $motif->binding_matrix->{'ic'}->[$i] . ",";
								}
								if ($i == $motif->length() - 1) {
									$infoCont .= $motif->binding_matrix->{'ic'}->[$i];
								}
							}
					
							unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]}{Difference} <= $totalScore) {
								$motifScoreHash{$motifName[0]}{Start} = $StartPos;
								$motifScoreHash{$motifName[0]}{Length} = $motiflength;
								$motifScoreHash{$motifName[0]}{DistToStart} = $DistToStart;
								$motifScoreHash{$motifName[0]}{SizeOfsubstitution} = $SizeOfsubstitution;
								$motifScoreHash{$motifName[0]}{Original} = $refScore;
								$motifScoreHash{$motifName[0]}{New} = max(@newBindingScores);
								$motifScoreHash{$motifName[0]}{Difference} = $totalScore;
								$motifScoreHash{$motifName[0]}{InfoContent} = $infoCont;
								$motifScoreHash{$motifName[0]}{ConsScore} = $consScore;
							}
						}
					}
	
					# For Insertions
					if (substr($vf->allele_string(), 0, 1) eq "-") {	
						if ($motif->strand() == 1) {
							# Calculating binding scores
							my $refString = $slice->seq();
							my $refScore = $motif->binding_matrix->relative_affinity($refString);
				
							# Creating a sequence including the insertion
							# A slice is made with 1000 bases in each direction of the motif
							my $insertionSlice = $slice_adaptor->fetch_by_region('chromosome', $vf->slice->seq_region_name(), $motif->seq_region_start() - 1000, $motif->seq_region_end() + 1000);
							my $insertionString = $insertionSlice->seq();
							my $insertion = substr($vf->allele_string(), 2);
							$insertionString = substr($insertionString, 0, 1000 + $vf->start() - $motif->seq_region_start()) . $insertion . substr($insertionString, 1000 + $vf->start() - $motif->seq_region_start());							
							my @newBindingScores = ();
							for (my $i = 0; $i < $motif->length + length($insertion) + 1; $i++) {
								my $altString = substr($insertionString, 1000 + $vf->start() - $motif->seq_region_start() - $motif->length() + $i, $motif->seq_region_end() - $motif->seq_region_start() + 1);
								push @newBindingScores, $motif->binding_matrix->relative_affinity($altString);
							}
							my $totalScore = max(@newBindingScores)-$refScore;
						
							# Motif start
							my $StartPos = $motif->seq_region_start();
					
							# Length of motif
							my $motiflength = $motif->binding_matrix->length();
					
							# Distance from variant to start of motif
							my $DistToStart = $vf->start() - $motif->seq_region_start();
					
							# Size of indel
							my @splitallele = split("/", $vf->allele_string());
							my $SizeOfIndel = length($splitallele[1]);
					
							# Name of motif
							my @motifName = $motif->display_label();
						
							# Conservation score calculation
							my $scores = $csscore_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $motiflength);
							my ($consScore) = calculateConservationScore($scores);
						
							# TFBS information content
							my $infoCont = ();
							for (my $i = 0; $i < $motif->length(); $i++) {
								if ($i != $motif->length() - 1) {
									$infoCont .= $motif->binding_matrix->{'ic'}->[$i] . ",";
								}
								if ($i == $motif->length() - 1) {
									$infoCont .= $motif->binding_matrix->{'ic'}->[$i];
								}
							}
					
							unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]}{Difference} <= $totalScore) {
								$motifScoreHash{$motifName[0]}{Start} = $StartPos;
								$motifScoreHash{$motifName[0]}{Length} = $motiflength;
								$motifScoreHash{$motifName[0]}{DistToStart} = $DistToStart;
								$motifScoreHash{$motifName[0]}{SizeOfIndel} = $SizeOfIndel;
								$motifScoreHash{$motifName[0]}{Original} = $refScore;
								$motifScoreHash{$motifName[0]}{New} = max(@newBindingScores);
								$motifScoreHash{$motifName[0]}{Difference} = $totalScore;
								$motifScoreHash{$motifName[0]}{InfoContent} = $infoCont;
								$motifScoreHash{$motifName[0]}{ConsScore} = $consScore;
							}
						}
			
						if ($motif->strand == -1) {
							# Calculating binding scores
							my $refString = PatternModule::RegAnno::createComplementarySequence($slice->seq);
							my $refScore = $motif->binding_matrix->relative_affinity($refString);
				
							# Creating a sequence including the insertion
							# A slice is made with 1000 bases in each direction of the motif
							my $insertionSlice = $slice_adaptor->fetch_by_region('chromosome', $vf->slice->seq_region_name(), $motif->seq_region_start() - 1000, $motif->seq_region_end() + 1000);
							my $insertionString = $insertionSlice->seq();
							my $insertion = substr($vf->allele_string(), 2);
							$insertionString = substr($insertionString, 0, 1000 + $vf->start() - $motif->seq_region_start()) . $insertion . substr($insertionString, 1000 + $vf->start() - $motif->seq_region_start());
							my @newBindingScores = ();
							for (my $i = 0; $i < $motif->length + length($insertion) + 1; $i++) {
								my $altString = PatternModule::RegAnno::createComplementarySequence(substr($insertionString, 1000 + $vf->start() - $motif->seq_region_start() - $motif->length() + $i, $motif->seq_region_end() - $motif->seq_region_start() + 1));
								push @newBindingScores, $motif->binding_matrix->relative_affinity($altString);
							}
							my $totalScore = max(@newBindingScores)-$refScore;
						
							# Motif start
							my $StartPos = $motif->seq_region_start();
						
							# Length of motif
							my $motiflength = $motif->binding_matrix->length();
					
							# Distance from variant to start of motif
							my $DistToStart = $motif->seq_region_end() - $vf->start();
					
							# Size of indel
							my @splitallele = split("/", $vf->allele_string());
							my $SizeOfIndel = length($splitallele[1]);
					
							# Name of motif
							my @motifName = $motif->display_label();
						
							# Conservation score calculation
							my $scores = $csscore_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $motiflength);
							my ($consScore) = calculateConservationScore($scores);
						
							# TFBS information content
							my $infoCont = ();
							for (my $i = 0; $i < $motif->length(); $i++) {
								if ($i != $motif->length() - 1) {
									$infoCont .= $motif->binding_matrix->{'ic'}->[$i] . ",";
								}
								if ($i == $motif->length() - 1) {
									$infoCont .= $motif->binding_matrix->{'ic'}->[$i];
								}
							}
					
							unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]}{Difference} <= $totalScore) {
								$motifScoreHash{$motifName[0]}{Start} = $StartPos;
								$motifScoreHash{$motifName[0]}{Length} = $motiflength;
								$motifScoreHash{$motifName[0]}{DistToStart} = $DistToStart;
								$motifScoreHash{$motifName[0]}{SizeOfIndel} = $SizeOfIndel;
								$motifScoreHash{$motifName[0]}{Original} = $refScore;
								$motifScoreHash{$motifName[0]}{New} = max(@newBindingScores);
								$motifScoreHash{$motifName[0]}{Difference} = $totalScore;
								$motifScoreHash{$motifName[0]}{InfoContent} = $infoCont;
								$motifScoreHash{$motifName[0]}{ConsScore} = $consScore;
							}
						}
					}
	
					# For Deletions
					if (substr($vf->allele_string(), -1) eq "-") {
						if ($motif->strand() == 1) {
					
							# Calculating binding scores
							my $refString = $slice->seq();
							my $refScore = $motif->binding_matrix->relative_affinity($refString);
			
							# Creating a sequence including the Deletion	
							# A slice is made with 1000 bases in each direction of the motif
							my $deletionSlice = $slice_adaptor->fetch_by_region('chromosome', $vf->slice->seq_region_name(), $motif->seq_region_start() - 1000, $motif->seq_region_end() + 1000);
							my $deletionString = $deletionSlice->seq();
							substr($deletionString, $vf->start() - $motif->seq_region_start() + 1000, $vf->end() - $vf->start() + 1) = "";
							my @newBindingScores = ();						
							for (my $i = 0; $i < $motif->length + 1; $i++) {
								my $altString = substr($deletionString, 1000 + $vf->start() - $motif->seq_region_start() - $motif->length() + $i, $motif->length());
								push @newBindingScores, $motif->binding_matrix->relative_affinity($altString);
							}
							my $totalScore = max(@newBindingScores)-$refScore;
						
							# Motif start
							my $StartPos = $motif->seq_region_start();
					
							# Length of motif
							my $motiflength = $motif->binding_matrix->length();
		
							# Distance from variant to start of motif
							my $DistToStart = -1;
							if ($vf->start() - $motif->seq_region_start() > -1) {
								$DistToStart = $vf->start() - $motif->seq_region_start();
							}
					
							# Size of indel
							my @splitallele = split("/", $vf->allele_string());
							my $SizeOfIndel = length($splitallele[0]);
					
							# If the deletion is of the entire transcription factor binding motif
							if ($vf->start() - $motif->seq_region_start() <= 0 && $motif->seq_region_end() - $vf->end() <= 0) {
								$SizeOfIndel = $motif->binding_matrix->length();
							}
					
							# If the deletion start/end both lie within the transcription factor binding motif
							if ($vf->start() - $motif->seq_region_start() >= 0 && $motif->seq_region_end() - $vf->end() >= 0) {
								$SizeOfIndel = length($splitallele[0]);
							}
					
							# If the deletion start lies before the motif start, and the deletion end lies within the motif
							if ($vf->start() - $motif->seq_region_start() <= 0 && $motif->seq_region_end() - $vf->end() >= 0) {
								$SizeOfIndel = $vf->end() - $motif->seq_region_start() + 1;
							}
					
							# If the deletion start lies within motif, and the deletion end lies outside motif
							if ($vf->start() - $motif->seq_region_start() >= 0 && $motif->seq_region_end() - $vf->end() <= 0) {
								$SizeOfIndel = $motif->seq_region_end() - $vf->start() + 1;
							}								
					
							# Name of motif
							my @motifName = $motif->display_label();
						
							# Conservation score calculation
							my $scores = $csscore_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $motiflength);
							my ($consScore) = calculateConservationScore($scores);
						
							# TFBS information content
							my $infoCont = ();
							for (my $i = 0; $i < $motif->length(); $i++) {
								if ($i != $motif->length() - 1) {
									$infoCont .= $motif->binding_matrix->{'ic'}->[$i] . ",";
								}
								if ($i == $motif->length() - 1) {
									$infoCont .= $motif->binding_matrix->{'ic'}->[$i];
								}
							}
					
							unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]}{Difference} <= $totalScore) {
								$motifScoreHash{$motifName[0]}{Start} = $StartPos;
								$motifScoreHash{$motifName[0]}{Length} = $motiflength;
								$motifScoreHash{$motifName[0]}{DistToStart} = $DistToStart;
								$motifScoreHash{$motifName[0]}{SizeOfIndel} = $SizeOfIndel;
								$motifScoreHash{$motifName[0]}{Original} = $refScore;
								$motifScoreHash{$motifName[0]}{New} = max(@newBindingScores);
								$motifScoreHash{$motifName[0]}{Difference} = $totalScore;
								$motifScoreHash{$motifName[0]}{InfoContent} = $infoCont;
								$motifScoreHash{$motifName[0]}{ConsScore} = $consScore;
							}
						}
						if ($motif->strand == -1) {
					
							# Calculating binding scores
							my $refString = PatternModule::RegAnno::createComplementarySequence($slice->seq());
							my $refScore = $motif->binding_matrix->relative_affinity($refString);
			
							# Creating a sequence including the Deletion
							# A slice is made with 1000 bases in each direction of the motif
							my $deletionSlice = $slice_adaptor->fetch_by_region('chromosome', $vf->slice->seq_region_name(), $motif->seq_region_start() - 1000, $motif->seq_region_end() + 1000);
							my $deletionString = $deletionSlice->seq();
							substr($deletionString, $vf->start - $motif->seq_region_start() + 1000, $vf->end() - $vf->start() + 1) = "";
							my @newBindingScores = ();						
							for (my $i = 0; $i < $motif->length() + 1; $i++) {
								my $altString = PatternModule::RegAnno::createComplementarySequence(substr($deletionString, 1000 + $vf->start() - $motif->seq_region_start() - $motif->length() + $i, $motif->length()));
								push @newBindingScores, $motif->binding_matrix->relative_affinity($altString);
							}
							my $totalScore = max(@newBindingScores)-$refScore;
						
							# Motif start
							my $StartPos = $motif->seq_region_start();
						
							# Length of motif
							my $motiflength = $motif->binding_matrix->length();

							# Distance from variant to start of motif
							my $DistToStart = -1;
							if ($motif->seq_region_end() - $vf->end() > -1) {
								$DistToStart = $motif->seq_region_end() - $vf->end();
							}
					
							# Size of indel
							my @splitallele = split("/", $vf->allele_string());
							my $SizeOfIndel = length($splitallele[0]);
					
							# If the deletion is of the entire transcription factor binding motif
							if ($vf->start() - $motif->seq_region_start() <= 0 && $motif->seq_region_end() - $vf->end() <= 0) {
								$SizeOfIndel = $motif->binding_matrix->length();
							}
					
							# If the deletion start/end both lie within the transcription factor binding motif
							if ($vf->start() - $motif->seq_region_start() >= 0 && $motif->seq_region_end() - $vf->end() >= 0) {
								$SizeOfIndel = length($splitallele[0]);
							}
					
							# If the deletion start lies before the motif start, and the deletion end lies within the motif
							if ($vf->start() - $motif->seq_region_start() <= 0 && $motif->seq_region_end() - $vf->end() >= 0) {
								$SizeOfIndel = $vf->end() - $motif->seq_region_start() + 1;
							}
					
							# If the deletion start lies within motif, and the deletion end lies outside motif
							if ($vf->start() - $motif->seq_region_start() >= 0 && $motif->seq_region_end() - $vf->end() <= 0) {
								$SizeOfIndel = $motif->seq_region_end() - $vf->start() + 1;
							}		
					
							# Name of motif
							my @motifName = $motif->display_label();
						
							# Conservation score calculation
							my $scores = $csscore_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $motiflength);
							my ($consScore) = calculateConservationScore($scores);
						
							# TFBS information content
							my $infoCont = ();
							for (my $i = 0; $i < $motif->length(); $i++) {
								if ($i != $motif->length() - 1) {
									$infoCont .= $motif->binding_matrix->{'ic'}->[$i] . ",";
								}
								if ($i == $motif->length() - 1) {
									$infoCont .= $motif->binding_matrix->{'ic'}->[$i];
								}
							}
					
							unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]}{Difference} <= $totalScore) {
								$motifScoreHash{$motifName[0]}{Start} = $StartPos;
								$motifScoreHash{$motifName[0]}{Length} = $motiflength;
								$motifScoreHash{$motifName[0]}{DistToStart} = $DistToStart;
								$motifScoreHash{$motifName[0]}{SizeOfIndel} = $SizeOfIndel;
								$motifScoreHash{$motifName[0]}{Original} = $refScore;
								$motifScoreHash{$motifName[0]}{New} = max(@newBindingScores);
								$motifScoreHash{$motifName[0]}{Difference} = $totalScore;
								$motifScoreHash{$motifName[0]}{InfoContent} = $infoCont;
								$motifScoreHash{$motifName[0]}{ConsScore} = $consScore;
							}
						}
					}
				}
			}
	
			# TFBS names are extracted
			@keys = sort keys %motifScoreHash;

			# A new simple hash with only TFBS name, and Difference between old and new binding scores are present
			for my $key (@keys) {
				$motifScoreSortingHash{$key} = $motifScoreHash{$key}{Difference};
			}
	
			# Keys are sorted such that the lowest Difference-score is first
			my @DifferenceSortedKeys = sort {$motifScoreSortingHash{$a} <=> $motifScoreSortingHash{$b}} keys %motifScoreSortingHash;
	
			# The TFBS information is added to the annotation string
			for my $Rkey (@DifferenceSortedKeys) {
				# Splitting TFBS name and matrix into two
				my @TFBSnames = split(":", $Rkey);
			
				$annotationString .= $TFBSnames[0] . "\t" . $TFBSnames[-1] . "\t" . $motifScoreHash{$Rkey}{Start} . "\t" . $motifScoreHash{$Rkey}{Length} . "\t" . $motifScoreHash{$Rkey}{DistToStart} . "\t" . $motifScoreHash{$Rkey}{SizeOfIndel} . "\t" . $motifScoreHash{$Rkey}{Original} . "\t" . $motifScoreHash{$Rkey}{New} . "\t" . $motifScoreHash{$Rkey}{InfoContent} . "\t" . $motifScoreHash{$Rkey}{ConsScore}  . "\t";
			}
	
			# The empty spaces in the annotation string are filled out with NA
			if (scalar @keys != 0) {
				for (my $i = 0; $i < 10 - scalar @keys; $i++) {
					$annotationString .= "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t";
				}
			}
		}	
	
	# If the variant overlaps with no TFBS, NA is filled in into the annotation string
	if (scalar @keys == 0) {
		$annotationString .= "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" ."NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t". "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t". "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t";
	}
	
	return $annotationString;
}

=head2 fetchMultiallelicVariantTFBSinfo

  Args       : "Multiallelic VariationFeature", "MotifFeature Adaptor", "Slice Adaptor"
  Example    : my $string = fetchVariantTFBSinfo($mvf, $motifFeature_adaptor, $slice_adaptor)
  Description: Creates a string with information on what transcription factor binding sites 
  			   overlap a variant, the conservation score of the TFBS and how the predicted 
  			   binding is affected by the variant
  Returntype : String

=cut

sub fetchMultiallelicVariantTFBSinfo {				
	my ($vf, $vf_list, $motifFeature_adaptor, $slice_adaptor, $csscore_adaptor, $mlss) = @_;	
	my @variationFeatureList = @$vf_list;
	my @allelicString = split("/", $vf->allele_string());
	my $variantSlice = $slice_adaptor->fetch_by_region('chromosome',$vf->slice->seq_region_name(),$vf->start(),$vf->start()-1+length($allelicString[0]));
	my @motif_features = @{$motifFeature_adaptor->fetch_all_by_Slice($variantSlice)};
	my @keys = ();
	
	my $annotationString = ();
	
	# Getting information on transcription factor binding sites that are hit by variants
	if (scalar @motif_features > 0) {
		my %motifScoreHash = ();
		my %motifScoreSortingHash = ();
	
		for my $varfeat (@variationFeatureList) {
			for my $motif (@motif_features) {
				if ($varfeat->start() <= $motif->seq_region_end() && $motif->seq_region_start() <= $varfeat->end()) {
					my $slice = $slice_adaptor->fetch_by_region('chromosome', $varfeat->slice->seq_region_name(), $motif->seq_region_start(), $motif->seq_region_end());				
									
					# For SNPs
					if ($vf->allele_string() !~ /-/ && length($vf->allele_string()) == 3) {
						if ($motif->strand() == 1) {
					
							# Calculating binding scores
							my $refString = $slice->seq();
							my $refScore = $motif->binding_matrix->relative_affinity($refString);
							my $snpIndex = $varfeat->start() - $slice->start();
							my $altString = $refString;
							substr($altString, $snpIndex, 1) = substr($varfeat->allele_string(), 2, 1);
							my $altScore = $motif->binding_matrix->relative_affinity($altString);
							my $totalScore = $altScore-$refScore;
							
							# Motif start
							my $StartPos = $motif->seq_region_start();
							
							# Length of motif
							my $motiflength = $motif->binding_matrix->length();
					
							# Distance from variant to start of motif
							my $DistToStart = $varfeat->start() - $motif->seq_region_start();
					
							# Size of indel
							my $SizeOfIndel = "NA";
				
							# Name of motif
							my @motifName = $motif->display_label();
							
							# Conservation score calculation
							my $scores = $csscore_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $motiflength);
							my ($consScore) = calculateConservationScore($scores);
							
							# TFBS information content
							my $infoCont = ();
							for (my $i = 0; $i < $motif->length(); $i++) {
								if ($i != $motif->length() - 1) {
									$infoCont .= $motif->binding_matrix->{'ic'}->[$i] . ",";
								}
								if ($i == $motif->length() - 1) {
									$infoCont .= $motif->binding_matrix->{'ic'}->[$i];
								}
							}
					
							# Placing information into Hash. If the same motif has two scores for a variant, the most severe score is stored
							unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]}{Difference} <= $totalScore) {
								$motifScoreHash{$motifName[0]}{Start} = $StartPos;
								$motifScoreHash{$motifName[0]}{Length} = $motiflength;
								$motifScoreHash{$motifName[0]}{DistToStart} = $DistToStart;
								$motifScoreHash{$motifName[0]}{SizeOfIndel} = $SizeOfIndel;
								$motifScoreHash{$motifName[0]}{Original} = $refScore;
								$motifScoreHash{$motifName[0]}{New} = $altScore;
								$motifScoreHash{$motifName[0]}{Difference} = $totalScore;
								$motifScoreHash{$motifName[0]}{InfoContent} = $infoCont;
								$motifScoreHash{$motifName[0]}{ConsScore} = $consScore;
							}
						}
			
						if ($motif->strand() == -1) {
					
							# Calculating binding scores
							my $refString = PatternModule::RegAnno::createComplementarySequence($slice->seq());
							my $refScore = $motif->binding_matrix->relative_affinity($refString);
							my $snpIndex = $varfeat->start() - $slice->start();
							my $altString = $slice->seq();
							substr($altString, $snpIndex, 1) = substr($varfeat->allele_string(), 2, 1);
							$altString = PatternModule::RegAnno::createComplementarySequence($altString);
							my $altScore = $motif->binding_matrix->relative_affinity($altString);
							my $totalScore = $altScore-$refScore;
							
							# Motif start
							my $StartPos = $motif->seq_region_start();
							
							# Length of motif
							my $motiflength = $motif->binding_matrix->length();
				
							# Distance from variant to start of motif
							my $DistToStart = $motif->seq_region_end() - $varfeat->start();
					
							# Size of indel
							my $SizeOfIndel = "NA";
					
							# Name of motif
							my @motifName = $motif->display_label();
							
							# Conservation score calculation
							my $scores = $csscore_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $motiflength);
							my ($consScore) = calculateConservationScore($scores);
							
							# TFBS information content
							my $infoCont = ();
							for (my $i = 0; $i < $motif->length(); $i++) {
								if ($i != $motif->length() - 1) {
									$infoCont .= $motif->binding_matrix->{'ic'}->[$i] . ",";
								}
								if ($i == $motif->length() - 1) {
									$infoCont .= $motif->binding_matrix->{'ic'}->[$i];
								}
							}
					
							unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]}{Difference} <= $totalScore) {
								$motifScoreHash{$motifName[0]}{Start} = $StartPos;
								$motifScoreHash{$motifName[0]}{Length} = $motiflength;
								$motifScoreHash{$motifName[0]}{DistToStart} = $DistToStart;
								$motifScoreHash{$motifName[0]}{SizeOfIndel} = $SizeOfIndel;
								$motifScoreHash{$motifName[0]}{Original} = $refScore;
								$motifScoreHash{$motifName[0]}{New} = $altScore;
								$motifScoreHash{$motifName[0]}{Difference} = $totalScore;
								$motifScoreHash{$motifName[0]}{InfoContent} = $infoCont;
								$motifScoreHash{$motifName[0]}{ConsScore} = $consScore;
							}
						}	
					}
	
				# For Substitutions
				if ($varfeat->allele_string() !~ /-/ && length($varfeat->allele_string()) > 3) {
					if ($motif->strand() == 1) {
					
						# Calculating binding scores
						my $refString = $slice->seq();
						my $refScore = $motif->binding_matrix->relative_affinity($refString);
			
						# Creating a sequence including the substitution	
						# A slice is made with 1000 bases in each direction of the motif
						my $deletionSlice = $slice_adaptor->fetch_by_region('chromosome', $varfeat->slice->seq_region_name(), $motif->seq_region_start() - 1000, $motif->seq_region_end() + 1000);
						my $deletionString = $deletionSlice->seq();
						my @varrefalt = split("/", $varfeat->allele_string());
						substr($deletionString, $varfeat->start() - $motif->seq_region_start() + 1000, $varfeat->end() - $varfeat->start() + 1) = $varrefalt[1];
						my @newBindingScores = ();						
						for (my $i = 0; $i < $motif->length + 1; $i++) {
							my $altString = substr($deletionString, 1000 + $varfeat->start() - $motif->seq_region_start() - $motif->length() + $i, $motif->length());
							push @newBindingScores, $motif->binding_matrix->relative_affinity($altString);
						}
						my $totalScore = max(@newBindingScores)-$refScore;
						
						# Motif start
						my $StartPos = $motif->seq_region_start();
						
						# Length of motif
						my $motiflength = $motif->binding_matrix->length();
		
						# Distance from variant to start of motif
						my $DistToStart = -1;
						if ($varfeat->start() - $motif->seq_region_start() > -1) {
							$DistToStart = $varfeat->start() - $motif->seq_region_start();
						}
					
						# Size of substitution
						my @splitallele = split("/", $varfeat->allele_string());
						my $SizeOfsubstitution = length($splitallele[0]);
					
						# If the substitution is of the entire transcription factor binding motif
						if ($varfeat->start() - $motif->seq_region_start() <= 0 && $motif->seq_region_end() - $varfeat->end() <= 0) {
							$SizeOfsubstitution = $motif->binding_matrix->length();
						}
					
						# If the substitution start/end both lie within the transcription factor binding motif
						if ($varfeat->start() - $motif->seq_region_start() >= 0 && $motif->seq_region_end() - $varfeat->end() >= 0) {
							$SizeOfsubstitution = length($splitallele[0]);
						}
					
						# If the substitution start lies before the motif start, and the substitution end lies within the motif
						if ($varfeat->start() - $motif->seq_region_start() <= 0 && $motif->seq_region_end() - $varfeat->end() >= 0) {
							$SizeOfsubstitution = $varfeat->end() - $motif->seq_region_start() + 1;
						}
					
						# If the substitution start lies within motif, and the substitution end lies outside motif
						if ($varfeat->start() - $motif->seq_region_start() >= 0 && $motif->seq_region_end() - $varfeat->end() <= 0) {
							$SizeOfsubstitution = $motif->seq_region_end() - $varfeat->start() + 1;
						}								
					
						# Name of motif
						my @motifName = $motif->display_label();
						
						# Conservation score calculation
						my $scores = $csscore_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $motiflength);
						my ($consScore) = calculateConservationScore($scores);
					
						# TFBS information content
						my $infoCont = ();
						for (my $i = 0; $i < $motif->length(); $i++) {
							if ($i != $motif->length() - 1) {
								$infoCont .= $motif->binding_matrix->{'ic'}->[$i] . ",";
							}
							if ($i == $motif->length() - 1) {
								$infoCont .= $motif->binding_matrix->{'ic'}->[$i];
							}
						}
					
						unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]}{Difference} <= $totalScore) {
							$motifScoreHash{$motifName[0]}{Start} = $StartPos;
							$motifScoreHash{$motifName[0]}{Length} = $motiflength;
							$motifScoreHash{$motifName[0]}{DistToStart} = $DistToStart;
							$motifScoreHash{$motifName[0]}{SizeOfsubstitution} = $SizeOfsubstitution;
							$motifScoreHash{$motifName[0]}{Original} = $refScore;
							$motifScoreHash{$motifName[0]}{New} = max(@newBindingScores);
							$motifScoreHash{$motifName[0]}{Difference} = $totalScore;
							$motifScoreHash{$motifName[0]}{InfoContent} = $infoCont;
							$motifScoreHash{$motifName[0]}{ConsScore} = $consScore;
						}
					}
					if ($motif->strand == -1) {
					
						# Calculating binding scores
						my $refString = PatternModule::RegAnno::createComplementarySequence($slice->seq());
						my $refScore = $motif->binding_matrix->relative_affinity($refString);
			
						# Creating a sequence including the substitution
						# A slice is made with 1000 bases in each direction of the motif
						my $deletionSlice = $slice_adaptor->fetch_by_region('chromosome', $varfeat->slice->seq_region_name(), $motif->seq_region_start() - 1000, $motif->seq_region_end() + 1000);
						my $deletionString = $deletionSlice->seq();
						my @varrefalt = split("/", $varfeat->allele_string());
						substr($deletionString, $varfeat->start - $motif->seq_region_start() + 1000, $varfeat->end() - $varfeat->start() + 1) = $varrefalt[1];
						my @newBindingScores = ();						
						for (my $i = 0; $i < $motif->length() + 1; $i++) {
							my $altString = PatternModule::RegAnno::createComplementarySequence(substr($deletionString, 1000 + $varfeat->start() - $motif->seq_region_start() - $motif->length() + $i, $motif->length()));
							push @newBindingScores, $motif->binding_matrix->relative_affinity($altString);
						}
						my $totalScore = max(@newBindingScores)-$refScore;
						
						# Motif start
						my $StartPos = $motif->seq_region_start();
						
						# Length of motif
						my $motiflength = $motif->binding_matrix->length();

						# Distance from variant to start of motif
						my $DistToStart = -1;
						if ($motif->seq_region_end() - $varfeat->end() > -1) {
							$DistToStart = $motif->seq_region_end() - $varfeat->end();
						}
					
						# Size of substitution
						my @splitallele = split("/", $varfeat->allele_string());
						my $SizeOfsubstitution = length($splitallele[0]);
					
						# If the substitution is of the entire transcription factor binding motif
						if ($varfeat->start() - $motif->seq_region_start() <= 0 && $motif->seq_region_end() - $varfeat->end() <= 0) {
							$SizeOfsubstitution = $motif->binding_matrix->length();
						}
					
						# If the substitution start/end both lie within the transcription factor binding motif
						if ($varfeat->start() - $motif->seq_region_start() >= 0 && $motif->seq_region_end() - $varfeat->end() >= 0) {
							$SizeOfsubstitution = length($splitallele[0]);
						}
					
						# If the substitution start lies before the motif start, and the substitution end lies within the motif
						if ($varfeat->start() - $motif->seq_region_start() <= 0 && $motif->seq_region_end() - $varfeat->end() >= 0) {
							$SizeOfsubstitution = $varfeat->end() - $motif->seq_region_start() + 1;
						}
					
						# If the substitution start lies within motif, and the substitution end lies outside motif
						if ($varfeat->start() - $motif->seq_region_start() >= 0 && $motif->seq_region_end() - $varfeat->end() <= 0) {
							$SizeOfsubstitution = $motif->seq_region_end() - $varfeat->start() + 1;
						}		
					
						# Name of motif
						my @motifName = $motif->display_label();
						
						# Conservation score calculation
						my $scores = $csscore_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $motiflength);
						my ($consScore) = calculateConservationScore($scores);
						
						# TFBS information content
						my $infoCont = ();
						for (my $i = 0; $i < $motif->length(); $i++) {
							if ($i != $motif->length() - 1) {
								$infoCont .= $motif->binding_matrix->{'ic'}->[$i] . ",";
							}
							if ($i == $motif->length() - 1) {
								$infoCont .= $motif->binding_matrix->{'ic'}->[$i];
							}
						}
					
						unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]}{Difference} <= $totalScore) {
							$motifScoreHash{$motifName[0]}{Start} = $StartPos;
							$motifScoreHash{$motifName[0]}{Length} = $motiflength;
							$motifScoreHash{$motifName[0]}{DistToStart} = $DistToStart;
							$motifScoreHash{$motifName[0]}{SizeOfsubstitution} = $SizeOfsubstitution;
							$motifScoreHash{$motifName[0]}{Original} = $refScore;
							$motifScoreHash{$motifName[0]}{New} = max(@newBindingScores);
							$motifScoreHash{$motifName[0]}{Difference} = $totalScore;
							$motifScoreHash{$motifName[0]}{InfoContent} = $infoCont;
							$motifScoreHash{$motifName[0]}{ConsScore} = $consScore;
						}
					}
				}				
	
					# For Insertions
					if (substr($varfeat->allele_string(), 0, 1) eq "-") {											
						if ($motif->strand() == 1) {
					
							# Calculating binding scores
							my $refString = $slice->seq();
							my $refScore = $motif->binding_matrix->relative_affinity($refString);
				
							# Creating a sequence including the insertion
							# A slice is made with 1000 bases in each direction of the motif
							my $insertionSlice = $slice_adaptor->fetch_by_region('chromosome', $varfeat->slice->seq_region_name(), $motif->seq_region_start() - 1000, $motif->seq_region_end() + 1000);
							my $insertionString = $insertionSlice->seq();
							my $insertion = substr($varfeat->allele_string(), 2);
							$insertionString = substr($insertionString, 0, 1000 + $varfeat->start() - $motif->seq_region_start()) . $insertion . substr($insertionString, 1000 + $varfeat->start() - $motif->seq_region_start());							
							my @newBindingScores = ();
							for (my $i = 0; $i < $motif->length + length($insertion) + 1; $i++) {
								my $altString = substr($insertionString, 1000 + $varfeat->start() - $motif->seq_region_start() - $motif->length() + $i, $motif->seq_region_end() - $motif->seq_region_start() + 1);
								push @newBindingScores, $motif->binding_matrix->relative_affinity($altString);
							}
							my $totalScore = max(@newBindingScores)-$refScore;
							
							# Motif start
							my $StartPos = $motif->seq_region_start();
							
							# Length of motif
							my $motiflength = $motif->binding_matrix->length();
					
							# Distance from variant to start of motif
							my $DistToStart = $varfeat->start() - $motif->seq_region_start();
					
							# Size of indel
							my @splitallele = split("/", $varfeat->allele_string());
							my $SizeOfIndel = length($splitallele[1]);
					
							# Name of motif
							my @motifName = $motif->display_label();
							
							# Conservation score calculation
							my $scores = $csscore_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $motiflength);
							my ($consScore) = calculateConservationScore($scores);
							
							# TFBS information content
							my $infoCont = ();
							for (my $i = 0; $i < $motif->length(); $i++) {
								if ($i != $motif->length() - 1) {
									$infoCont .= $motif->binding_matrix->{'ic'}->[$i] . ",";
								}
								if ($i == $motif->length() - 1) {
									$infoCont .= $motif->binding_matrix->{'ic'}->[$i];
								}
							}
					
							unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]}{Difference} <= $totalScore) {
								$motifScoreHash{$motifName[0]}{Start} = $StartPos;
								$motifScoreHash{$motifName[0]}{Length} = $motiflength;
								$motifScoreHash{$motifName[0]}{DistToStart} = $DistToStart;
								$motifScoreHash{$motifName[0]}{SizeOfIndel} = $SizeOfIndel;
								$motifScoreHash{$motifName[0]}{Original} = $refScore;
								$motifScoreHash{$motifName[0]}{New} = max(@newBindingScores);
								$motifScoreHash{$motifName[0]}{Difference} = $totalScore;
								$motifScoreHash{$motifName[0]}{InfoContent} = $infoCont;
								$motifScoreHash{$motifName[0]}{ConsScore} = $consScore;
							}
						}
			
						if ($motif->strand == -1) {
					
							# Calculating binding scores
							my $refString = PatternModule::RegAnno::createComplementarySequence($slice->seq);
							my $refScore = $motif->binding_matrix->relative_affinity($refString);
				
							# Creating a sequence including the insertion
							# A slice is made with 1000 bases in each direction of the motif
							my $insertionSlice = $slice_adaptor->fetch_by_region('chromosome', $varfeat->slice->seq_region_name(), $motif->seq_region_start() - 1000, $motif->seq_region_end() + 1000);
							my $insertionString = $insertionSlice->seq();
							my $insertion = substr($varfeat->allele_string(), 2);
							$insertionString = substr($insertionString, 0, 1000 + $varfeat->start() - $motif->seq_region_start()) . $insertion . substr($insertionString, 1000 + $varfeat->start() - $motif->seq_region_start());
							my @newBindingScores = ();
							for (my $i = 0; $i < $motif->length + length($insertion) + 1; $i++) {
								my $altString = PatternModule::RegAnno::createComplementarySequence(substr($insertionString, 1000 + $varfeat->start() - $motif->seq_region_start() - $motif->length() + $i, $motif->seq_region_end() - $motif->seq_region_start() + 1));
								push @newBindingScores, $motif->binding_matrix->relative_affinity($altString);
							}
							my $totalScore = max(@newBindingScores)-$refScore;
							
							# Motif start
							my $StartPos = $motif->seq_region_start();
							
							# Length of motif
							my $motiflength = $motif->binding_matrix->length();
					
							# Distance from variant to start of motif
							my $DistToStart = $motif->seq_region_end() - $varfeat->start();
					
							# Size of indel
							my @splitallele = split("/", $varfeat->allele_string());
							my $SizeOfIndel = length($splitallele[1]);
					
							# Name of motif
							my @motifName = $motif->display_label();
							
							# Conservation score calculation
							my $scores = $csscore_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $motiflength);
							my ($consScore) = calculateConservationScore($scores);
							
							# TFBS information content
							my $infoCont = ();
							for (my $i = 0; $i < $motif->length(); $i++) {
								if ($i != $motif->length() - 1) {
									$infoCont .= $motif->binding_matrix->{'ic'}->[$i] . ",";
								}
								if ($i == $motif->length() - 1) {
									$infoCont .= $motif->binding_matrix->{'ic'}->[$i];
								}
							}
					
							unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]}{Difference} <= $totalScore) {
								$motifScoreHash{$motifName[0]}{Start} = $StartPos;
								$motifScoreHash{$motifName[0]}{Length} = $motiflength;
								$motifScoreHash{$motifName[0]}{DistToStart} = $DistToStart;
								$motifScoreHash{$motifName[0]}{SizeOfIndel} = $SizeOfIndel;
								$motifScoreHash{$motifName[0]}{Original} = $refScore;
								$motifScoreHash{$motifName[0]}{New} = max(@newBindingScores);
								$motifScoreHash{$motifName[0]}{Difference} = $totalScore;
								$motifScoreHash{$motifName[0]}{InfoContent} = $infoCont;
								$motifScoreHash{$motifName[0]}{ConsScore} = $consScore;
							}
						}
					}
	
					# For Deletions
					if (substr($varfeat->allele_string(), -1) eq "-") {
						if ($motif->strand() == 1) {
					
							# Calculating binding scores
							my $refString = $slice->seq();
							my $refScore = $motif->binding_matrix->relative_affinity($refString);
			
							# Creating a sequence including the Deletion	
							# A slice is made with 1000 bases in each direction of the motif
							my $deletionSlice = $slice_adaptor->fetch_by_region('chromosome', $varfeat->slice->seq_region_name(), $motif->seq_region_start() - 1000, $motif->seq_region_end() + 1000);
							my $deletionString = $deletionSlice->seq();
							substr($deletionString, $varfeat->start() - $motif->seq_region_start() + 1000, $varfeat->end() - $varfeat->start() + 1) = "";
							my @newBindingScores = ();						
							for (my $i = 0; $i < $motif->length + 1; $i++) {
								my $altString = substr($deletionString, 1000 + $varfeat->start() - $motif->seq_region_start() - $motif->length() + $i, $motif->length());
								push @newBindingScores, $motif->binding_matrix->relative_affinity($altString);
							}
							my $totalScore = max(@newBindingScores)-$refScore;
							
							# Motif start
							my $StartPos = $motif->seq_region_start();
							
							# Length of motif
							my $motiflength = $motif->binding_matrix->length();
		
							# Distance from variant to start of motif
							my $DistToStart = -1;
							if ($varfeat->start() - $motif->seq_region_start() > -1) {
								$DistToStart = $varfeat->start() - $motif->seq_region_start();
							}
					
							# Size of indel
							my @splitallele = split("/", $varfeat->allele_string());
							my $SizeOfIndel = length($splitallele[0]);
					
							# If the deletion is of the entire transcription factor binding motif
							if ($varfeat->start() - $motif->seq_region_start() <= 0 && $motif->seq_region_end() - $varfeat->end() <= 0) {
								$SizeOfIndel = $motif->binding_matrix->length();
							}
					
							# If the deletion start/end both lie within the transcription factor binding motif
							if ($varfeat->start() - $motif->seq_region_start() >= 0 && $motif->seq_region_end() - $varfeat->end() >= 0) {
								$SizeOfIndel = length($splitallele[0]);
							}
					
							# If the deletion start lies before the motif start, and the deletion end lies within the motif
							if ($varfeat->start() - $motif->seq_region_start() <= 0 && $motif->seq_region_end() - $varfeat->end() >= 0) {
								$SizeOfIndel = $varfeat->end() - $motif->seq_region_start() + 1;
							}
					
							# If the deletion start lies within motif, and the deletion end lies outside motif
							if ($varfeat->start() - $motif->seq_region_start() >= 0 && $motif->seq_region_end() - $varfeat->end() <= 0) {
								$SizeOfIndel = $motif->seq_region_end() - $varfeat->start() + 1;
							}								
					
							# Name of motif
							my @motifName = $motif->display_label();
							
							# Conservation score calculation
							my $scores = $csscore_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $motiflength);
							my ($consScore) = calculateConservationScore($scores);
							
							# TFBS information content
							my $infoCont = ();
							for (my $i = 0; $i < $motif->length(); $i++) {
								if ($i != $motif->length() - 1) {
									$infoCont .= $motif->binding_matrix->{'ic'}->[$i] . ",";
								}
								if ($i == $motif->length() - 1) {
									$infoCont .= $motif->binding_matrix->{'ic'}->[$i];
								}
							}
					
							unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]}{Difference} <= $totalScore) {
								$motifScoreHash{$motifName[0]}{Start} = $StartPos;
								$motifScoreHash{$motifName[0]}{Length} = $motiflength;
								$motifScoreHash{$motifName[0]}{DistToStart} = $DistToStart;
								$motifScoreHash{$motifName[0]}{SizeOfIndel} = $SizeOfIndel;
								$motifScoreHash{$motifName[0]}{Original} = $refScore;
								$motifScoreHash{$motifName[0]}{New} = max(@newBindingScores);
								$motifScoreHash{$motifName[0]}{Difference} = $totalScore;
								$motifScoreHash{$motifName[0]}{InfoContent} = $infoCont;
								$motifScoreHash{$motifName[0]}{ConsScore} = $consScore;
							}
						}
						if ($motif->strand == -1) {
					
							# Calculating binding scores
							my $refString = PatternModule::RegAnno::createComplementarySequence($slice->seq());
							my $refScore = $motif->binding_matrix->relative_affinity($refString);
			
							# Creating a sequence including the Deletion
							# A slice is made with 1000 bases in each direction of the motif
							my $deletionSlice = $slice_adaptor->fetch_by_region('chromosome', $varfeat->slice->seq_region_name(), $motif->seq_region_start() - 1000, $motif->seq_region_end() + 1000);
							my $deletionString = $deletionSlice->seq();
							substr($deletionString, $varfeat->start - $motif->seq_region_start() + 1000, $varfeat->end() - $varfeat->start() + 1) = "";
							my @newBindingScores = ();						
							for (my $i = 0; $i < $motif->length() + 1; $i++) {
								my $altString = PatternModule::RegAnno::createComplementarySequence(substr($deletionString, 1000 + $varfeat->start() - $motif->seq_region_start() - $motif->length() + $i, $motif->length()));
								push @newBindingScores, $motif->binding_matrix->relative_affinity($altString);
							}
							my $totalScore = max(@newBindingScores)-$refScore;
							
							# Motif start
							my $StartPos = $motif->seq_region_start();
							
							# Length of motif
							my $motiflength = $motif->binding_matrix->length();

							# Distance from variant to start of motif
							my $DistToStart = -1;
							if ($motif->seq_region_end() - $varfeat->end() > -1) {
								$DistToStart = $motif->seq_region_end() - $varfeat->end();
							}
					
							# Size of indel
							my @splitallele = split("/", $varfeat->allele_string());
							my $SizeOfIndel = length($splitallele[0]);
					
							# If the deletion is of the entire transcription factor binding motif
							if ($varfeat->start() - $motif->seq_region_start() <= 0 && $motif->seq_region_end() - $varfeat->end() <= 0) {
								$SizeOfIndel = $motif->binding_matrix->length();
							}
					
							# If the deletion start/end both lie within the transcription factor binding motif
							if ($varfeat->start() - $motif->seq_region_start() >= 0 && $motif->seq_region_end() - $varfeat->end() >= 0) {
								$SizeOfIndel = length($splitallele[0]);
							}
					
							# If the deletion start lies before the motif start, and the deletion end lies within the motif
							if ($varfeat->start() - $motif->seq_region_start() <= 0 && $motif->seq_region_end() - $varfeat->end() >= 0) {
								$SizeOfIndel = $varfeat->end() - $motif->seq_region_start() + 1;
							}
					
							# If the deletion start lies within motif, and the deletion end lies outside motif
							if ($varfeat->start() - $motif->seq_region_start() >= 0 && $motif->seq_region_end() - $varfeat->end() <= 0) {
								$SizeOfIndel = $motif->seq_region_end() - $varfeat->start() + 1;
							}		
					
							# Name of motif
							my @motifName = $motif->display_label();
							
							# Conservation score calculation
							my $scores = $csscore_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $motiflength);
							my ($consScore) = calculateConservationScore($scores);
							
							# TFBS information content
							my $infoCont = ();
							for (my $i = 0; $i < $motif->length(); $i++) {
								if ($i != $motif->length() - 1) {
									$infoCont .= $motif->binding_matrix->{'ic'}->[$i] . ",";
								}
								if ($i == $motif->length() - 1) {
									$infoCont .= $motif->binding_matrix->{'ic'}->[$i];
								}
							}
					
							unless (exists $motifScoreHash{$motifName[0]} and $motifScoreHash{$motifName[0]}{Difference} <= $totalScore) {
								$motifScoreHash{$motifName[0]}{Start} = $StartPos;
								$motifScoreHash{$motifName[0]}{Length} = $motiflength;
								$motifScoreHash{$motifName[0]}{DistToStart} = $DistToStart;
								$motifScoreHash{$motifName[0]}{SizeOfIndel} = $SizeOfIndel;
								$motifScoreHash{$motifName[0]}{Original} = $refScore;
								$motifScoreHash{$motifName[0]}{New} = max(@newBindingScores);
								$motifScoreHash{$motifName[0]}{Difference} = $totalScore;
								$motifScoreHash{$motifName[0]}{InfoContent} = $infoCont;
								$motifScoreHash{$motifName[0]}{ConsScore} = $consScore;
							}
						}
					}
				}
			}
		}				
	
		# TFBS names are extracted
		@keys = sort keys %motifScoreHash;
	
		# A new simple hash with only TFBS name, and Difference between old and new binding scores are present
		for my $key (@keys) {
			$motifScoreSortingHash{$key} = $motifScoreHash{$key}{Difference};
		}
	
		# Keys are sorted such that the lowest Difference-score is first
		my @DifferenceSortedKeys = sort {$motifScoreSortingHash{$a} <=> $motifScoreSortingHash{$b}} keys %motifScoreSortingHash;
	
		# The TFBS information is added to the annotation string
		for my $Rkey (@DifferenceSortedKeys) {
			# Splitting TFBS name and matrix into two
			my @TFBSnames = split(":", $Rkey);
			
			$annotationString .= $TFBSnames[0] . "\t" . $TFBSnames[-1] . "\t" . $motifScoreHash{$Rkey}{Start} . "\t" . $motifScoreHash{$Rkey}{Length} . "\t" . $motifScoreHash{$Rkey}{DistToStart} . "\t" . $motifScoreHash{$Rkey}{SizeOfIndel} . "\t" . $motifScoreHash{$Rkey}{Original} . "\t" . $motifScoreHash{$Rkey}{New} . "\t" . $motifScoreHash{$Rkey}{InfoContent} . "\t" . $motifScoreHash{$Rkey}{ConsScore}  . "\t";
		}
	
		# The empty spaces in the annotation string are filled out with NA
		if (scalar @keys != 0) {
			for (my $i = 0; $i < 10 - scalar @keys; $i++) {
				$annotationString .= "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t";
			}
		}
	}	

	# If the variant overlaps with no TFBS, NA is filled in into the annotation string
	if (scalar @keys == 0) {
		$annotationString .= "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" ."NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t". "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t". "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t" . "NA" . "\t";
	}
			
	return $annotationString;
}

=head2 calculateConservationScore

  Args       : "Conservation Scores"
  Example    : my $consScore = calculateConservationScore($scores);
  Description: Calculates the average conservation score given a combination of conservation
  			   scores and missing values
  Returntype : Integer

=cut

sub calculateConservationScore {
	my ($scores) = @_;
	my $sum = 0;
	my $count = 0;
	my $diff_score;
	foreach my $score (@$scores) {
		if (defined $score->diff_score) {
			$sum += $score->diff_score;
			$count++;
		}
	}	
	if ($count > 0){
		$diff_score = $sum / $count;
	}
	if (defined $diff_score){
		return $diff_score;
	}
	else {
		return "NA";
	}
}

=head2 stringToHashConverter

  Args       : "Tab-separated String"
  Example    : my %TFBSHash = stringToHashConverter($TFBSstring, \%variantInformation);
  Description: Places information from a tab-separated string into a hash on the appropriate
  			   places
  Returntype : Hash

=cut

sub stringToHashConverter {
	my ($string, $hash) = @_;
	my %infoHash = %{$hash};
	my @infoArray = split("\t", $string);
	my $counter = 1;
	
	for (my $i = 0; $i < 100; $i = $i + 10) {
		# Creating variable key names for the ten different TFBS annotations
		my $TFBS1string = "ENSTFBSANNO_TFBS" . $counter . "_NAME";
		my $TFBS2string = "ENSTFBSANNO_TFBS" . $counter . "_MATRIX";
		my $TFBS3string = "ENSTFBSANNO_TFBS" . $counter . "_POS";
		my $TFBS4string = "ENSTFBSANNO_TFBS" . $counter . "_LENGTH";
		my $TFBS5string = "ENSTFBSANNO_TFBS" . $counter . "_DIST_TO_START";
		my $TFBS6string = "ENSTFBSANNO_TFBS" . $counter . "_SIZE_OF_INDEL";
		my $TFBS7string = "ENSTFBSANNO_TFBS" . $counter . "_REF_SCORE";
		my $TFBS8string = "ENSTFBSANNO_TFBS" . $counter . "_ALT_SCORE";
		my $TFBS9string = "ENSTFBSANNO_TFBS" . $counter . "_INFO_CONTENT";
		my $TFBS10string = "ENSTFBSANNO_TFBS" . $counter . "_CONS_SCORE";
		
		# Placing information in array into hash		
		if ($infoArray[$i] ne "NA") {
			$infoHash{$TFBS1string} = $infoArray[$i];
		}
		if ($infoArray[$i+1] ne "NA") {
			$infoHash{$TFBS2string} = $infoArray[$i+1];
		}
		if ($infoArray[$i+2] ne "NA") {
			$infoHash{$TFBS3string} = $infoArray[$i+2];
		}
		if ($infoArray[$i+3] ne "NA") {
			$infoHash{$TFBS4string} = $infoArray[$i+3];
		}
		if ($infoArray[$i+4] ne "NA") {
			$infoHash{$TFBS5string} = $infoArray[$i+4];
		}
		if ($infoArray[$i+5] ne "NA") {
			$infoHash{$TFBS6string} = $infoArray[$i+5];
		}
		if ($infoArray[$i+6] ne "NA") {
			$infoHash{$TFBS7string} = $infoArray[$i+6];
		}
		if ($infoArray[$i+7] ne "NA") {
			$infoHash{$TFBS8string} = $infoArray[$i+7];
		}
		if ($infoArray[$i+8] ne "NA") {
			$infoHash{$TFBS9string} = $infoArray[$i+8];
		}
		if ($infoArray[$i+9] ne "NA") {
			$infoHash{$TFBS10string} = $infoArray[$i+9];
		}
		
		# Counts what TFBS we are at
		$counter = $counter + 1;
	}

	return \%infoHash;
}

return 1;