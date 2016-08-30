package PatternModule::PredictAnno;

use strict;
use warnings;

=head2 createVariantPredictHash

  Args       : None
  Example    : my %variantInformation = createVariantPredictHash()
  Description: Creates a hash that stores information on the variant's predicted PolyPhen
  			   and SIFT consequences
  Returntype : Hash

=cut

sub createVariantPredictHash {
	my %data = ('ENSPREDICTANNO_POLYPHEN_TRANSCRIPT'			=> "NA",
				'ENSPREDICTANNO_POLYPHEN_TRANSCRIPT_NAME'		=> "NA",
				'ENSPREDICTANNO_POLYPHEN_PREDICTION'			=> "NA",
				'ENSPREDICTANNO_POLYPHEN_SCORE'					=> ".",
				'ENSPREDICTANNO_SIFT_TRANSCRIPT'				=> "NA",
				'ENSPREDICTANNO_SIFT_TRANSCRIPT_NAME'			=> "NA",
				'ENSPREDICTANNO_SIFT_PREDICTION'				=> "NA",
				'ENSPREDICTANNO_SIFT_SCORE'						=> ".",
	);	
	return %data;
}

=head2 predictConsequences

  Args       : "VariationFeature"
  Example    : my ($SIFT_transcript, $SIFT_prediction, $SIFT_score, $PolyPhen_transcript, 
  			   $PolyPhen_prediction, $PolyPhen_score) = 
  			   PatternModule::PredictAnno::predictConsequences($vf);
  Description: For each variant the most severe PolyPhen and SIFT predictions are annotated
  	 		   with the transcript the predictions are made on. For transcripts with 
  	 		   identical scores, the longest transcript is annotated
  Returntype : List

=cut

sub predictConsequences {
	my ($vf) = @_;
	
	# Ranking of SIFT predictions for sorting consequences
	my %SIFT_prediction = (	'deleterious'		=> "1",
							'tolerated'			=> "2",
			);	
	
	# Ranking of PolyPhen predictions for sorting consequences
	my %PolyPhen_prediction = (	'probably_damaging'		=> "1",
								'possibly_damaging'		=> "2",
								'benign'				=> "3",
								'unknown'				=> "4",
	);	
	
	# Creating a hash for storing and sorting scores, such that only the most severe score is kept
	my %scoreHash = createVariantPredictHash();
	
	# All transcriptvariations are tested for scores on all transcripts
	my @TranscriptVariations = @{$vf->get_all_TranscriptVariations()};
	for my $tv (@TranscriptVariations) {
		
		# The alleles are scored
		my @TranscriptVariationAlleles = @{$tv->get_all_TranscriptVariationAlleles()};
		for my $tva (@TranscriptVariationAlleles) {
			########
			# SIFT #
			########
			
			# Calculating SIFT predictions
			my $SIFTScore = $tva->sift_score();
			my $SIFTPredict = $tva->sift_prediction();
			
			# If SIFT Score exists and there are no lower scores in the hash, update hash
			if ($SIFTScore) {
				
				# If there is not entry in the hash yet, the SIFT data is placed there
				if ($scoreHash{"ENSPREDICTANNO_SIFT_TRANSCRIPT"} eq "NA") {
					$scoreHash{"ENSPREDICTANNO_SIFT_TRANSCRIPT"} = $tv;
					$scoreHash{"ENSPREDICTANNO_SIFT_TRANSCRIPT_NAME"} = $tv->transcript->stable_id();
					$scoreHash{"ENSPREDICTANNO_SIFT_PREDICTION"} = $SIFTPredict;
					$scoreHash{"ENSPREDICTANNO_SIFT_SCORE"} = $SIFTScore;
				}
				
				# If has already contains SIFT data, the transcript will only be added if 
				# it has a more severe score or prediction
				if ($scoreHash{"ENSPREDICTANNO_SIFT_TRANSCRIPT"} ne "NA") {
					if ($scoreHash{"ENSPREDICTANNO_SIFT_SCORE"} > $SIFTScore or $SIFT_prediction{$scoreHash{"ENSPREDICTANNO_SIFT_PREDICTION"}} > $SIFT_prediction{$SIFTPredict}) {
						$scoreHash{"ENSPREDICTANNO_SIFT_TRANSCRIPT"} = $tv;
						$scoreHash{"ENSPREDICTANNO_SIFT_TRANSCRIPT_NAME"} = $tv->transcript->stable_id();
						$scoreHash{"ENSPREDICTANNO_SIFT_PREDICTION"} = $SIFTPredict;
						$scoreHash{"ENSPREDICTANNO_SIFT_SCORE"} = $SIFTScore;
					}
					
					# If both have same prediction and score the longest transcript is chosen
					if ($scoreHash{"ENSPREDICTANNO_SIFT_SCORE"} == $SIFTScore and $SIFT_prediction{$scoreHash{"ENSPREDICTANNO_SIFT_PREDICTION"}} == $SIFT_prediction{$SIFTPredict}) {
						if ($scoreHash{"ENSPREDICTANNO_SIFT_TRANSCRIPT"}->transcript->length() < $tv->transcript->length()) {
							$scoreHash{"ENSPREDICTANNO_SIFT_TRANSCRIPT"} = $tv;
							$scoreHash{"ENSPREDICTANNO_SIFT_TRANSCRIPT_NAME"} = $tv->transcript->stable_id();
							$scoreHash{"ENSPREDICTANNO_SIFT_PREDICTION"} = $SIFTPredict;
							$scoreHash{"ENSPREDICTANNO_SIFT_SCORE"} = $SIFTScore;
						}
					}
				}				
			}
			
			############
			# PolyPhen #
			############
			
			# Calculating PolyPhen predictions
			my $PolyPhenScore = $tva->polyphen_score();
			my $PolyPhenPredict = $tva->polyphen_prediction();
			
			# Removing white space from annotation
			if ($PolyPhenPredict) {
				$PolyPhenPredict =~ tr/ /_/;
			}
			
			# If SIFT Score exists and there are no lower scores in the hash, update hash
			if ($PolyPhenScore) {
				
				# If there is not entry in the hash yet, the POLYPHEN data is placed there
				if ($scoreHash{"ENSPREDICTANNO_POLYPHEN_TRANSCRIPT"} eq "NA") {
					$scoreHash{"ENSPREDICTANNO_POLYPHEN_TRANSCRIPT"} = $tv;
					$scoreHash{"ENSPREDICTANNO_POLYPHEN_TRANSCRIPT_NAME"} = $tv->transcript->stable_id();
					$scoreHash{"ENSPREDICTANNO_POLYPHEN_PREDICTION"} = $PolyPhenPredict;
					$scoreHash{"ENSPREDICTANNO_POLYPHEN_SCORE"} = $PolyPhenScore;
				}
				
				# If has already contains POLYPHEN data, the transcript will only be added if 
				# it has a more severe score or prediction
				if ($scoreHash{"ENSPREDICTANNO_POLYPHEN_TRANSCRIPT"} ne "NA") {
					if ($scoreHash{"ENSPREDICTANNO_POLYPHEN_SCORE"} < $PolyPhenScore or $PolyPhen_prediction{$scoreHash{"ENSPREDICTANNO_POLYPHEN_PREDICTION"}} > $PolyPhen_prediction{$PolyPhenPredict}) {
						$scoreHash{"ENSPREDICTANNO_POLYPHEN_TRANSCRIPT"} = $tv;
						$scoreHash{"ENSPREDICTANNO_POLYPHEN_TRANSCRIPT_NAME"} = $tv->transcript->stable_id();
						$scoreHash{"ENSPREDICTANNO_POLYPHEN_PREDICTION"} = $PolyPhenPredict;
						$scoreHash{"ENSPREDICTANNO_POLYPHEN_SCORE"} = $PolyPhenScore;
					}
					
					# If both have same prediction and score the longest transcript is chosen
					if ($scoreHash{"ENSPREDICTANNO_POLYPHEN_SCORE"} == $PolyPhenScore and $PolyPhen_prediction{$scoreHash{"ENSPREDICTANNO_POLYPHEN_PREDICTION"}} == $PolyPhen_prediction{$PolyPhenPredict}) {
						if ($scoreHash{"ENSPREDICTANNO_POLYPHEN_TRANSCRIPT"}->transcript->length() < $tv->transcript->length()) {
							$scoreHash{"ENSPREDICTANNO_POLYPHEN_TRANSCRIPT"} = $tv;
							$scoreHash{"ENSPREDICTANNO_POLYPHEN_TRANSCRIPT_NAME"} = $tv->transcript->stable_id();
							$scoreHash{"ENSPREDICTANNO_POLYPHEN_PREDICTION"} = $PolyPhenPredict;
							$scoreHash{"ENSPREDICTANNO_POLYPHEN_SCORE"} = $PolyPhenScore;
						}
					}
				}
			}		
		}
	}
	return %scoreHash;
}

return 1;