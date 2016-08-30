#!/usr/bin/perl -w

# Modules and packages
use strict;
use warnings;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::ApiVersion;
use PatternModule::GeneAnno;
use PatternModule::RegAnno;
use PatternModule::CleaningVariants;
use PatternModule::TFBSAnno;
use PatternModule::PredictAnno;
use List::Util qw(min);
use List::Util qw(max);
use List::MoreUtils qw(firstidx);

#########################
##### Documentation #####
#########################

=head1 NAME

Variant Effect Annotation Tool

=cut

=head1 DESCRIPTION

This script is designed to annotate the effects and positions of variants in coding, 
non-coding and regulatory regions. If a variant does not fall within a specific region, 
the nearest corresponding region is found.

=cut

=head1 SYNOPSIS

=head3 USAGE

perl PatternAnnotator.pl -in vcf_file_input.vcf -out vcf_file_output.vcf -all

=head3 INPUT

The script takes as input a VCF file that has been processed to split the variants that 
contain multiple alternative alleles into multiple lines, with one alternative allele per 
line for those individuals that are heterozygous for the reference. For variants that are
heterozygous for two alternative alleles, the lines remain multi allelic (see details in 
the associated script "MultiAllelicSplit.pl"). When the alternative alleles are "complex" 
(i.e. SNPs, INS or DEL in the formats AA/AG, AT/ATT, ATT/AT), the variants are further 
"cleaned" by a module in the script, in order to correct for sequence overlaps and other 
format issues (See CleaningVariants.pm).

=head3 ARGUMENTS
 
-all														
By choosing "-all" you get annotations of protein coding, long non-coding, short 
non-coding and regulatory regions at the same time. If you only need a subset of the 
annotations the "-pc", "-lnc", "-snc", -tfbs and "-reg" options can be used either individually
or together.

-pc														
By choosing "-pc" your annotations will only be of protein coding genes. The genes and
transcripts that are chosen for your variants are determined in the following way: 1. The 
transcript with the "Most Severe Variant Consequence" is determined using ENSEMBLs 
prioritised consequence table (see consequence table in GeneAnno Module) 2. If several 
transcripts have the same consequence the canonical transcript is chosen 3. If the 
canonical transcript is not among the transcripts, the longest transcript is chosen. 4. 
From the most severely hit transcript the gene is determined. The position of the 
variant within the transcript/gene is annotated together with information on the length of
the gene and transcript. If the variant does not overlap with a gene, the nearest gene 
(genes within 5000 bases) is annotated. If you want to expand the interval within which 
closeby genes are found, you can used the option "-gene_interval". 

-lnc														
By choosing "-lnc" your annotations will only be of long non-coding genes. The genes and
transcripts that are chosen for your variants are determined in the same way as for 
protein coding genes/transcripts, and the same option is available. Additionally the 
transcript biotype is annotated.

-snc														
By choosing "-snc" your annotations will only be of short non-coding genes. The genes and
transcripts that are chosen for your variants are determined in the same way as for 
protein coding genes/transcripts, and the same option is available. Additionally the 
transcript biotype is annotated.

-reg														
By choosing "-reg" your annotations will only be of regulatory regions. We are annotating 
the transcription factor binding motifs separately from the other regulatory features 
(regions defined by ENSEMBL and referred to as RegulatoryFeature). Regulatory regions 
differ from cell type to cell type, which is why the annotations have been divided into 
three different groups of cells. MultiCell (see ENSEMBL's definition), ImmunoCell (Cell 
lines in the immune system - see RegAnno module for details) and NeuroCell (Cell lines in 
the central nervous system - see RegAnno module for details). The following information is 
annotated: 1. Whether a variant is located within a RegulatoryFeature (A region defined 
from multi cell evidence such as DNaseI hypersensitivity, transcription factor and 
polymerase enrichments (chipSeq) and histone modifications) 3. Position of variant within 
the RegulatoryFeature. 4. For each group of cell lines the evidence for the RegulatoryFeature 
is annotated 5. If a transcription factor binding motif is hit by a variant, the motif and 
the effect on the binding of the new motif (relative to that of the original reference 
sequence) is annotated. If multiple motifs for the same transcription factor are affected 
by the same variant, the most severely affected motif is annotated 6. If no RegulatoryFeature 
or transcription factor binding motif is affected, both the closest regulatory element and 
the closest transcription factor binding motif within 5000 bases is found (this can be 
adjusted using the options "-reg_interval" and "-tfbs_interval").

-tfbs
By choosing "-tfbs" your annotations will only be of transcription factor binding sites. 
The TFBS/TFBSs that the variant overlaps are annotated with their position weight matrices
with information content, their position relative to the start of the TFBS, the how great
a part of a deletion has removed of the TFBS, the gerp conservation score of the TFBS, and
predicted binding scores of the reference sequence and the alternative sequence. For 
multiallelic variants the most severe effect on binding is annotated. Annotations can handle
overlaps of up to ten TFBSs.

-gene_interval														
You can change the default value (5000 bases upstream and downstream of variant) for finding 
the nearest gene by using this option in the following way "-gene_interval 10000", where the 
number can be any number.

-reg_interval														
You can change the default value (5000 bases upstream and downstream of variant) for finding 
the nearest RegulatoryFeature by using this option in the following way "-reg_interval 
10000", where the number can be any number.

-tfbs_interval														
You can change the default value (5000 bases upstream and downstream of variant) for finding 
the nearest transcription factor binding motif by using this option in the following way 
"-tfbs_interval 10000", where the number can be any number.

-fix														
In fix mode only the lines containing the annotation IS_MULTI_ALLELIC=YES are re-annotated.
This option can be used for re-annotating the multi allelic variants after the original 
vcf-files have been split by individual, thus ensuring that the most severe consequences 
in a specific individual are annotated. This option can only be used in combination with
either of the -all, -pc, -lnc, -snc, -reg, -tfbs options.

=head3 RETURN		
												
The output is an annotatated vcf-file. The script fixes the format of those variants present
in the "complex" format (i.e. resulting from multisample calling) as explained previously. 
See detailed documentation on the cleaning methods (CleaningVariants.pm).

=head3 OTHER		
												
The script requires that the ENSEMBL API is installed on your system, as it is from this
database the annotations are pulled. For details on the installation see ENSEMBL's 
official website

=head1 FEEDBACK

Esben Eickhardt esbeneickhardt@biomed.au.dk

=head1 AUTHORS

Esben Eickhardt - Email E<lt>esbeneickhardt@biomed.au.dkE<gt>

Francesco Lescai - Email E<lt>lescai@biomed.au.dkE<gt>

=cut

##############################################
##### Connecting to the ENSEMBL database #####
##############################################
# For using online databases
# my $registry = 'Bio::EnsEMBL::Registry';
#
# $registry->load_registry_from_db(
#     -host => 'ensembldb.ensembl.org',
#     -user => 'anonymous'
# );

# For using local databases
my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
   -host => 'vm7',
   -user => 'anonymous',
   -database => 'human',
   -db_version => '75'
);
die ("Can't initialise registry") if (!$registry);

# Getting adaptors
my $cs_adaptor = $registry->get_adaptor( 'Human', 'Core', 'CoordSystem' );
my $slice_adaptor = $registry->get_adaptor( 'human', 'core', 'slice');
my $variant_feature_adaptor = $registry->get_adaptor('human', 'variation', 'variationfeature');
my $transcript_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Transcript' );
my $gene_adaptor  = $registry->get_adaptor( 'Human', 'Core', 'Gene' );
my $regfeat_adaptor = $registry->get_adaptor('Human', 'Funcgen', 'RegulatoryFeature');
my $motifFeature_adaptor = $registry->get_adaptor('Human', 'funcgen', 'motiffeature');
my $csscore_adaptor = $registry->get_adaptor("Multi", 'compara', 'ConservationScore');
my $mlss_adaptor = $registry->get_adaptor("Multi", "compara", "MethodLinkSpeciesSet");

# Get method_link_species_set object for gerp conservation scores for mammals
my $mlss = $mlss_adaptor->fetch_by_method_link_type_species_set_name("GERP_CONSERVATION_SCORE", "mammals");
throw("Unable to find method_link_species_set") if (!defined($mlss));

###############################################
##### Infolines to be added to vcf-header #####
###############################################

# For protein coding genes
my $info_line1 =  '##INFO=<ID=ENSPCANNO_VARIANT_OVERLAPS_WITH_GENE,Number=1,Type=String,Description="Ensembl API Boolean telling whether the variant falls within a protein coding gene">';
my $info_line2 =  '##INFO=<ID=ENSPCANNO_GENE,Number=1,Type=String,Description="Ensembl API Gene with most severely hit protein coding transcript">';
my $info_line3 =  '##INFO=<ID=ENSPCANNO_GENE_LENGTH,Number=1,Type=Integer,Description="Ensembl API Length of the gene with most severely hit protein coding transcript">';
my $info_line4 =  '##INFO=<ID=ENSPCANNO_GENE_NUMBER_OF_TRANSCRIPTS,Number=1,Type=Integer,Description="Ensembl API The number of transcripts the gene with most severely hit protein coding transcript has">';
my $info_line5 =  '##INFO=<ID=ENSPCANNO_GENE_DISTANCE_TO_START,Number=1,Type=Integer,Description="Ensembl API Distance from variant to gene start site">';
my $info_line6 =  '##INFO=<ID=ENSPCANNO_GENE_DISTANCE_TO_END,Number=1,Type=Integer,Description="Ensembl API Distance from variant to gene end site">';
my $info_line7 =  '##INFO=<ID=ENSPCANNO_TRANSCRIPT,Number=1,Type=String,Description="Ensembl API Most severely hit protein coding transcript decided by ENSEMBLs ranking of consequences, whether the transcript is canonical and whether it is the longest of the transcripts">';
my $info_line8 =  '##INFO=<ID=ENSPCANNO_TRANSCRIPT_BIOTYPE,Number=1,Type=String,Description="Ensembl API Transcript biotype is protein coding for all genes overlapping with a gene, for upstream variants various biotypes are possible">';
my $info_line9 =  '##INFO=<ID=ENSPCANNO_TRANSCRIPT_LENGTH,Number=1,Type=String,Description="Ensembl API Length of the most severely hit protein coding transcript">';
my $info_line10 =  '##INFO=<ID=ENSPCANNO_TRANSCRIPT_DISTANCE_TO_START,Number=1,Type=Integer,Description="Ensembl API variant distance to transcript start site is calculated">';
my $info_line11 =  '##INFO=<ID=ENSPCANNO_TRANSCRIPT_DISTANCE_TO_END,Number=1,Type=Integer,Description="Ensembl API variant distance to transcript end site is calculated">';
my $info_line12 =  '##INFO=<ID=ENSPCANNO_TRANSCRIPT_CONSEQUENCE,Number=1,Type=String,Description="Ensembl API Consequence on most severely hit protein coding transcript">';

# For long non-coding genes
my $info_line13 =  '##INFO=<ID=ENSLNCANNO_GENE,Number=1,Type=String,Description="Ensembl API Gene with most severely hit long non-coding transcript">';
my $info_line14 =  '##INFO=<ID=ENSLNCANNO_GENE_LENGTH,Number=1,Type=Integer,Description="Ensembl API Length of the gene with most severely hit long non-coding transcript">';
my $info_line15 =  '##INFO=<ID=ENSLNCANNO_GENE_NUMBER_OF_TRANSCRIPTS,Number=1,Type=Integer,Description="Ensembl API The number of transcripts the gene with most severely hit long non-coding transcript has">';
my $info_line16 =  '##INFO=<ID=ENSLNCANNO_GENE_DISTANCE_TO_START,Number=1,Type=Integer,Description="Ensembl API Distance from variant to gene start site">';
my $info_line17 =  '##INFO=<ID=ENSLNCANNO_GENE_DISTANCE_TO_END,Number=1,Type=Integer,Description="Ensembl API Distance from variant to gene end site">';
my $info_line18 =  '##INFO=<ID=ENSLNCANNO_TRANSCRIPT,Number=1,Type=String,Description="Ensembl API Most severely hit long non-coding transcript decided by ENSEMBLs ranking of consequences, whether the transcript is canonical and whether it is the longest of the transcripts">';
my $info_line19 =  '##INFO=<ID=ENSLNCANNO_TRANSCRIPT_BIOTYPE,Number=1,Type=String,Description="Ensembl API Transcript biotype">';
my $info_line20 =  '##INFO=<ID=ENSLNCANNO_TRANSCRIPT_LENGTH,Number=1,Type=Integer,Description="Ensembl API Length of the most severely hit long non-coding transcript">';
my $info_line21 =  '##INFO=<ID=ENSLNCANNO_TRANSCRIPT_DISTANCE_TO_START,Number=1,Type=Integer,Description="Ensembl API If the variant falls within an exon the distance to transcript start site is calculated">';
my $info_line22 =  '##INFO=<ID=ENSLNCANNO_TRANSCRIPT_DISTANCE_TO_END,Number=1,Type=Integer,Description="Ensembl API If the variant falls within an exon the distance to transcript end site is calculated">';
my $info_line23 =  '##INFO=<ID=ENSLNCANNO_TRANSCRIPT_CONSEQUENCE,Number=1,Type=String,Description="Ensembl API Consequence on most severely hit long non-coding transcript">';

# For short non-coding genes
my $info_line24 =  '##INFO=<ID=ENSSNCANNO_GENE,Number=1,Type=String,Description="Ensembl API Gene with most severely hit short non-coding transcript">';
my $info_line25 =  '##INFO=<ID=ENSSNCANNO_GENE_LENGTH,Number=1,Type=Integer,Description="Ensembl API Length of the gene with most severely hit short non-coding transcript">';
my $info_line26 =  '##INFO=<ID=ENSSNCANNO_GENE_NUMBER_OF_TRANSCRIPTS,Number=1,Type=Integer,Description="Ensembl API The number of transcripts the gene with most severely hit short non-coding transcript has">';
my $info_line27 =  '##INFO=<ID=ENSSNCANNO_GENE_DISTANCE_TO_START,Number=1,Type=Integer,Description="Ensembl API Distance from variant to gene start site">';
my $info_line28 =  '##INFO=<ID=ENSSNCANNO_GENE_DISTANCE_TO_END,Number=1,Type=Integer,Description="Ensembl API Distance from variant to gene end site">';
my $info_line29 =  '##INFO=<ID=ENSSNCANNO_TRANSCRIPT,Number=1,Type=String,Description="Ensembl API Most severely hit short non-coding transcript decided by ENSEMBLs ranking of consequences, whether the transcript is canonical and whether it is the shortest of the transcripts">';
my $info_line30 =  '##INFO=<ID=ENSSNCANNO_TRANSCRIPT_BIOTYPE,Number=1,Type=String,Description="Ensembl API Transcript biotype">';
my $info_line31 =  '##INFO=<ID=ENSSNCANNO_TRANSCRIPT_LENGTH,Number=1,Type=Integer,Description="Ensembl API Length of the most severely hit short non-coding transcript">';
my $info_line32 =  '##INFO=<ID=ENSSNCANNO_TRANSCRIPT_DISTANCE_TO_START,Number=1,Type=Integer,Description="Ensembl API If the variant falls within an exon the distance to transcript start site is calculated">';
my $info_line33 =  '##INFO=<ID=ENSSNCANNO_TRANSCRIPT_DISTANCE_TO_END,Number=1,Type=Integer,Description="Ensembl API If the variant falls within an exon the distance to transcript end site is calculated">';
my $info_line34 =  '##INFO=<ID=ENSSNCANNO_TRANSCRIPT_CONSEQUENCE,Number=1,Type=String,Description="Ensembl API Consequence on most severely hit short non-coding transcript">';

# For regulatory regions
my $info_line35 =  '##INFO=<ID=ENSREGANNO_VARIANT_OVERLAPS_WITH_REGULATORY,Number=1,Type=String,Description="Boolean telling whether the variant falls within a multicell regulatory feature">';
my $info_line36 =  '##INFO=<ID=ENSREGANNO_REGULATORY_FEATURE,Number=1,Type=String,Description="Ensembl API ID of the regulatory feature, which is the same for all cell types">';
my $info_line37 =  '##INFO=<ID=ENSREGANNO_REGULATORY_FEATURE_TYPE,Number=1,Type=String,Description="Ensembl API Type of the regulatory feature for multicell (differs per cell type)">';
my $info_line38 =  '##INFO=<ID=ENSREGANNO_DISTANCE_TO_START,Number=1,Type=Integer,Description="Ensembl API Distance from variant to regulatory feature start">';
my $info_line39 =  '##INFO=<ID=ENSREGANNO_DISTANCE_TO_END,Number=1,Type=Integer,Description="Ensembl API Distance from variant to regulatory feature end">';
my $info_line40 =  '##INFO=<ID=ENSREGANNO_CONSEQUENCE,Number=1,Type=String,Description="Ensembl API consequence of variant of regulatory feature">';
my $info_line41 =  '##INFO=<ID=ENSREGANNO_MULTICELL_TFBS_HIT,Number=1,Type=String,Description="Ensembl API Is a motif hit by the variant in MULTICELL">';
my $info_line42 =  '##INFO=<ID=ENSREGANNO_MULTICELL_TFBS,Number=1,Type=String,Description="Ensembl API MultiCell TFBS hit by the variant (or nearest TFBS), if multiple TFBS are hit annotation is Multi">';
my $info_line43 =  '##INFO=<ID=ENSREGANNO_MULTICELL_MULTITFBS,Number=.,Type=String,Description="Ensembl API If multiple MultiCell TFBS are hit by variant, the TFBSs are listed here">';
my $info_line44 =  '##INFO=<ID=ENSREGANNO_MULTICELL_TFBS_SCORE,Number=1,Type=Float,Description="Ensembl API Score of affinity loss/gain for TFBS measured as affinity score of ref sequence divided by score of alt sequence">';
my $info_line45 =  '##INFO=<ID=ENSREGANNO_MULTICELL_MULTITFBS_SCORE,Number=.,Type=Float,Description="Ensembl API Scores of affinity loss/gain for MULTITFBS measured as affinity score of ref sequence divided my score of alt sequence">';
my $info_line46 =  '##INFO=<ID=ENSREGANNO_MULTICELL_EVIDENCE,Number=1,Type=String,Description="Ensembl API MultiCell Evidence for regulatory region, if there are several types of evidence annotation is multi">';
my $info_line47 =  '##INFO=<ID=ENSREGANNO_MULTICELL_MULTIEVIDENCE,Number=.,Type=String,Description="Ensembl API If multiple MultiCell Evidence, the Evidence is listed here">';
my $info_line48 =  '##INFO=<ID=ENSREGANNO_NEUROCELL_TFBS_HIT,Number=1,Type=String,Description="Ensembl API Is a motif hit by the variant in NEUROCELL">';
my $info_line49 =  '##INFO=<ID=ENSREGANNO_NEUROCELL_TFBS,Number=1,Type=String,Description="Ensembl API NeuroCell TFBS hit by the variant (or nearest TFBS), if multiple TFBS are hit annotation is Multi">';
my $info_line50 =  '##INFO=<ID=ENSREGANNO_NEUROCELL_MULTITFBS,Number=.,Type=String,Description="Ensembl API If multiple NeuroCell TFBS are hit by variant, the TFBSs are listed here">';
my $info_line51 =  '##INFO=<ID=ENSREGANNO_NEUROCELL_TFBS_SCORE,Number=1,Type=Float,Description="Ensembl API Score of affinity loss/gain for TFBS measured as affinity score of ref sequence divided my score of alt sequence">';
my $info_line52 =  '##INFO=<ID=ENSREGANNO_NEUROCELL_MULTITFBS_SCORE,Number=.,Type=Float,Description="Ensembl API Scores of affinity loss/gain for MULTITFBS measured as affinity score of ref sequence divided my score of alt sequence">';
my $info_line53 =  '##INFO=<ID=ENSREGANNO_NEUROCELL_EVIDENCE,Number=1,Type=String,Description="Ensembl API NeuroCell Evidence for regulatory region, if there are several types of evidence annotation is multi">';
my $info_line54 =  '##INFO=<ID=ENSREGANNO_NEUROCELL_MULTIEVIDENCE,Number=.,Type=String,Description="Ensembl API If multiple NeuroCell Evidence, the Evidence is listed here">';
my $info_line55 =  '##INFO=<ID=ENSREGANNO_IMMUNOCELL_TFBS_HIT,Number=1,Type=String,Description="Ensembl API Is a motif hit by the variant in IMMUNOCELL">';
my $info_line56 =  '##INFO=<ID=ENSREGANNO_IMMUNOCELL_TFBS,Number=1,Type=String,Description="Ensembl API ImmunoCell TFBS hit by the variant (or nearest TFBS), if multiple TFBS are hit annotation is Multi">';
my $info_line57 =  '##INFO=<ID=ENSREGANNO_IMMUNOCELL_MULTITFBS,Number=.,Type=String,Description="Ensembl API If multiple ImmunoCell TFBS are hit by variant, the TFBSs are listed here">';
my $info_line58 =  '##INFO=<ID=ENSREGANNO_IMMUNOCELL_TFBS_SCORE,Number=1,Type=Float,Description="Ensembl API Score of affinity loss/gain for TFBS measured as affinity score of ref sequence divided my score of alt sequence">';
my $info_line59 =  '##INFO=<ID=ENSREGANNO_IMMUNOCELL_MULTITFBS_SCORE,Number=.,Type=Float,Description="Ensembl API Scores of affinity loss/gain for MULTITFBS measured as affinity score of ref sequence divided my score of alt sequence">';
my $info_line60 =  '##INFO=<ID=ENSREGANNO_IMMUNOCELL_EVIDENCE,Number=1,Type=String,Description="Ensembl API ImmunoCell Evidence for regulatory region, if there are several types of evidence annotation is multi">';
my $info_line61 =  '##INFO=<ID=ENSREGANNO_IMMUNOCELL_MULTIEVIDENCE,Number=.,Type=String,Description="Ensembl API If multiple ImmunoCell Evidence, the Evidence is listed here">';

# For transcription factor binding sites
my $info_line62 =  '##INFO=<ID=ENSTFBSANNO_TFBS1_NAME,Number=1,Type=String,Description="Ensembl API Name of the overlapping TFBS">';
my $info_line63 =  '##INFO=<ID=ENSTFBSANNO_TFBS1_MATRIX,Number=1,Type=String,Description="Ensembl API Name of position weight matrix belonging to the TFBS">';
my $info_line64 =  '##INFO=<ID=ENSTFBSANNO_TFBS1_POS,Number=1,Type=Integer,Description="Ensembl API Position of variant within the TFBS">';
my $info_line65 =  '##INFO=<ID=ENSTFBSANNO_TFBS1_LENGTH,Number=1,Type=Integer,Description="Ensembl API Length of the TFBS">';
my $info_line66 =  '##INFO=<ID=ENSTFBSANNO_TFBS1_DIST_TO_START,Number=1,Type=Integer,Description="Ensembl API Distance from variant to start of TFBS">';
my $info_line67 =  '##INFO=<ID=ENSTFBSANNO_TFBS1_SIZE_OF_INDEL,Number=1,Type=Integer,Description="Ensembl API Size of the insertion or size of the deleted part of the TFBS">';
my $info_line68 =  '##INFO=<ID=ENSTFBSANNO_TFBS1_REF_SCORE,Number=1,Type=Float,Description="Ensembl API Score of affinity loss/gain for TFBS measured as affinity score of ref sequence relative to the most optimal binding">';
my $info_line69 =  '##INFO=<ID=ENSTFBSANNO_TFBS1_ALT_SCORE,Number=1,Type=Float,Description="Ensembl API Scores of affinity loss/gain for TFBS measured as affinity score of alt sequence relative to the most optimal binding">';
my $info_line70 =  '##INFO=<ID=ENSTFBSANNO_TFBS1_INFO_CONTENT,Number=1,Type=String,Description="Ensembl API information content across the TFBS">';
my $info_line71 =  '##INFO=<ID=ENSTFBSANNO_TFBS1_CONS_SCORE,Number=1,Type=Float,Description="Ensembl API GERP conservation score for the entire TFBS">';
my $info_line72 =  '##INFO=<ID=ENSTFBSANNO_TFBS2_NAME,Number=1,Type=String,Description="Ensembl API Name of the overlapping TFBS">';
my $info_line73 =  '##INFO=<ID=ENSTFBSANNO_TFBS2_MATRIX,Number=1,Type=String,Description="Ensembl API Name of position weight matrix belonging to the TFBS">';
my $info_line74 =  '##INFO=<ID=ENSTFBSANNO_TFBS2_POS,Number=1,Type=Integer,Description="Ensembl API Position of variant within the TFBS">';
my $info_line75 =  '##INFO=<ID=ENSTFBSANNO_TFBS2_LENGTH,Number=1,Type=Integer,Description="Ensembl API Length of the TFBS">';
my $info_line76 =  '##INFO=<ID=ENSTFBSANNO_TFBS2_DIST_TO_START,Number=1,Type=Integer,Description="Ensembl API Distance from variant to start of TFBS">';
my $info_line77 =  '##INFO=<ID=ENSTFBSANNO_TFBS2_SIZE_OF_INDEL,Number=1,Type=Integer,Description="Ensembl API Size of the insertion or size of the deleted part of the TFBS">';
my $info_line78 =  '##INFO=<ID=ENSTFBSANNO_TFBS2_REF_SCORE,Number=1,Type=Float,Description="Ensembl API Score of affinity loss/gain for TFBS measured as affinity score of ref sequence relative to the most optimal binding">';
my $info_line79 =  '##INFO=<ID=ENSTFBSANNO_TFBS2_ALT_SCORE,Number=1,Type=Float,Description="Ensembl API Scores of affinity loss/gain for TFBS measured as affinity score of alt sequence relative to the most optimal binding">';
my $info_line80 =  '##INFO=<ID=ENSTFBSANNO_TFBS2_INFO_CONTENT,Number=1,Type=String,Description="Ensembl API information content across the TFBS">';
my $info_line81 =  '##INFO=<ID=ENSTFBSANNO_TFBS2_CONS_SCORE,Number=1,Type=Float,Description="Ensembl API GERP conservation score for the entire TFBS">';
my $info_line82 =  '##INFO=<ID=ENSTFBSANNO_TFBS3_NAME,Number=1,Type=String,Description="Ensembl API Name of the overlapping TFBS">';
my $info_line83 =  '##INFO=<ID=ENSTFBSANNO_TFBS3_MATRIX,Number=1,Type=String,Description="Ensembl API Name of position weight matrix belonging to the TFBS">';
my $info_line84 =  '##INFO=<ID=ENSTFBSANNO_TFBS3_POS,Number=1,Type=Integer,Description="Ensembl API Position of variant within the TFBS">';
my $info_line85 =  '##INFO=<ID=ENSTFBSANNO_TFBS3_LENGTH,Number=1,Type=Integer,Description="Ensembl API Length of the TFBS">';
my $info_line86 =  '##INFO=<ID=ENSTFBSANNO_TFBS3_DIST_TO_START,Number=1,Type=Integer,Description="Ensembl API Distance from variant to start of TFBS">';
my $info_line87 =  '##INFO=<ID=ENSTFBSANNO_TFBS3_SIZE_OF_INDEL,Number=1,Type=Integer,Description="Ensembl API Size of the insertion or size of the deleted part of the TFBS">';
my $info_line88 =  '##INFO=<ID=ENSTFBSANNO_TFBS3_REF_SCORE,Number=1,Type=Float,Description="Ensembl API Score of affinity loss/gain for TFBS measured as affinity score of ref sequence relative to the most optimal binding">';
my $info_line89 =  '##INFO=<ID=ENSTFBSANNO_TFBS3_ALT_SCORE,Number=1,Type=Float,Description="Ensembl API Scores of affinity loss/gain for TFBS measured as affinity score of alt sequence relative to the most optimal binding">';
my $info_line90 =  '##INFO=<ID=ENSTFBSANNO_TFBS3_INFO_CONTENT,Number=1,Type=String,Description="Ensembl API information content across the TFBS">';
my $info_line91 =  '##INFO=<ID=ENSTFBSANNO_TFBS3_CONS_SCORE,Number=1,Type=Float,Description="Ensembl API GERP conservation score for the entire TFBS">';
my $info_line92 =  '##INFO=<ID=ENSTFBSANNO_TFBS4_NAME,Number=1,Type=String,Description="Ensembl API Name of the overlapping TFBS">';
my $info_line93 =  '##INFO=<ID=ENSTFBSANNO_TFBS4_MATRIX,Number=1,Type=String,Description="Ensembl API Name of position weight matrix belonging to the TFBS">';
my $info_line94 =  '##INFO=<ID=ENSTFBSANNO_TFBS4_POS,Number=1,Type=Integer,Description="Ensembl API Position of variant within the TFBS">';
my $info_line95 =  '##INFO=<ID=ENSTFBSANNO_TFBS4_LENGTH,Number=1,Type=Integer,Description="Ensembl API Length of the TFBS">';
my $info_line96 =  '##INFO=<ID=ENSTFBSANNO_TFBS4_DIST_TO_START,Number=1,Type=Integer,Description="Ensembl API Distance from variant to start of TFBS">';
my $info_line97 =  '##INFO=<ID=ENSTFBSANNO_TFBS4_SIZE_OF_INDEL,Number=1,Type=Integer,Description="Ensembl API Size of the insertion or size of the deleted part of the TFBS">';
my $info_line98 =  '##INFO=<ID=ENSTFBSANNO_TFBS4_REF_SCORE,Number=1,Type=Float,Description="Ensembl API Score of affinity loss/gain for TFBS measured as affinity score of ref sequence relative to the most optimal binding">';
my $info_line99 =  '##INFO=<ID=ENSTFBSANNO_TFBS4_ALT_SCORE,Number=1,Type=Float,Description="Ensembl API Scores of affinity loss/gain for TFBS measured as affinity score of alt sequence relative to the most optimal binding">';
my $info_line100 =  '##INFO=<ID=ENSTFBSANNO_TFBS4_INFO_CONTENT,Number=1,Type=String,Description="Ensembl API information content across the TFBS">';
my $info_line101 =  '##INFO=<ID=ENSTFBSANNO_TFBS4_CONS_SCORE,Number=1,Type=Float,Description="Ensembl API GERP conservation score for the entire TFBS">';
my $info_line102 =  '##INFO=<ID=ENSTFBSANNO_TFBS5_NAME,Number=1,Type=String,Description="Ensembl API Name of the overlapping TFBS">';
my $info_line103 =  '##INFO=<ID=ENSTFBSANNO_TFBS5_MATRIX,Number=1,Type=String,Description="Ensembl API Name of position weight matrix belonging to the TFBS">';
my $info_line104 =  '##INFO=<ID=ENSTFBSANNO_TFBS5_POS,Number=1,Type=Integer,Description="Ensembl API Position of variant within the TFBS">';
my $info_line105 =  '##INFO=<ID=ENSTFBSANNO_TFBS5_LENGTH,Number=1,Type=Integer,Description="Ensembl API Length of the TFBS">';
my $info_line106 =  '##INFO=<ID=ENSTFBSANNO_TFBS5_DIST_TO_START,Number=1,Type=Integer,Description="Ensembl API Distance from variant to start of TFBS">';
my $info_line107 =  '##INFO=<ID=ENSTFBSANNO_TFBS5_SIZE_OF_INDEL,Number=1,Type=Integer,Description="Ensembl API Size of the insertion or size of the deleted part of the TFBS">';
my $info_line108 =  '##INFO=<ID=ENSTFBSANNO_TFBS5_REF_SCORE,Number=1,Type=Float,Description="Ensembl API Score of affinity loss/gain for TFBS measured as affinity score of ref sequence relative to the most optimal binding">';
my $info_line109 =  '##INFO=<ID=ENSTFBSANNO_TFBS5_ALT_SCORE,Number=1,Type=Float,Description="Ensembl API Scores of affinity loss/gain for TFBS measured as affinity score of alt sequence relative to the most optimal binding">';
my $info_line110 =  '##INFO=<ID=ENSTFBSANNO_TFBS5_INFO_CONTENT,Number=1,Type=String,Description="Ensembl API information content across the TFBS">';
my $info_line111 =  '##INFO=<ID=ENSTFBSANNO_TFBS5_CONS_SCORE,Number=1,Type=Float,Description="Ensembl API GERP conservation score for the entire TFBS">';
my $info_line112 =  '##INFO=<ID=ENSTFBSANNO_TFBS6_NAME,Number=1,Type=String,Description="Ensembl API Name of the overlapping TFBS">';
my $info_line113 =  '##INFO=<ID=ENSTFBSANNO_TFBS6_MATRIX,Number=1,Type=String,Description="Ensembl API Name of position weight matrix belonging to the TFBS">';
my $info_line114 =  '##INFO=<ID=ENSTFBSANNO_TFBS6_POS,Number=1,Type=Integer,Description="Ensembl API Position of variant within the TFBS">';
my $info_line115 =  '##INFO=<ID=ENSTFBSANNO_TFBS6_LENGTH,Number=1,Type=Integer,Description="Ensembl API Length of the TFBS">';
my $info_line116 =  '##INFO=<ID=ENSTFBSANNO_TFBS6_DIST_TO_START,Number=1,Type=Integer,Description="Ensembl API Distance from variant to start of TFBS">';
my $info_line117 =  '##INFO=<ID=ENSTFBSANNO_TFBS6_SIZE_OF_INDEL,Number=1,Type=Integer,Description="Ensembl API Size of the insertion or size of the deleted part of the TFBS">';
my $info_line118 =  '##INFO=<ID=ENSTFBSANNO_TFBS6_REF_SCORE,Number=1,Type=Float,Description="Ensembl API Score of affinity loss/gain for TFBS measured as affinity score of ref sequence relative to the most optimal binding">';
my $info_line119 =  '##INFO=<ID=ENSTFBSANNO_TFBS6_ALT_SCORE,Number=1,Type=Float,Description="Ensembl API Scores of affinity loss/gain for TFBS measured as affinity score of alt sequence relative to the most optimal binding">';
my $info_line120 =  '##INFO=<ID=ENSTFBSANNO_TFBS6_INFO_CONTENT,Number=1,Type=String,Description="Ensembl API information content across the TFBS">';
my $info_line121 =  '##INFO=<ID=ENSTFBSANNO_TFBS6_CONS_SCORE,Number=1,Type=Float,Description="Ensembl API GERP conservation score for the entire TFBS">';
my $info_line122 =  '##INFO=<ID=ENSTFBSANNO_TFBS7_NAME,Number=1,Type=String,Description="Ensembl API Name of the overlapping TFBS">';
my $info_line123 =  '##INFO=<ID=ENSTFBSANNO_TFBS7_MATRIX,Number=1,Type=String,Description="Ensembl API Name of position weight matrix belonging to the TFBS">';
my $info_line124 =  '##INFO=<ID=ENSTFBSANNO_TFBS7_POS,Number=1,Type=Integer,Description="Ensembl API Position of variant within the TFBS">';
my $info_line125 =  '##INFO=<ID=ENSTFBSANNO_TFBS7_LENGTH,Number=1,Type=Integer,Description="Ensembl API Length of the TFBS">';
my $info_line126 =  '##INFO=<ID=ENSTFBSANNO_TFBS7_DIST_TO_START,Number=1,Type=Integer,Description="Ensembl API Distance from variant to start of TFBS">';
my $info_line127 =  '##INFO=<ID=ENSTFBSANNO_TFBS7_SIZE_OF_INDEL,Number=1,Type=Integer,Description="Ensembl API Size of the insertion or size of the deleted part of the TFBS">';
my $info_line128 =  '##INFO=<ID=ENSTFBSANNO_TFBS7_REF_SCORE,Number=1,Type=Float,Description="Ensembl API Score of affinity loss/gain for TFBS measured as affinity score of ref sequence relative to the most optimal binding">';
my $info_line129 =  '##INFO=<ID=ENSTFBSANNO_TFBS7_ALT_SCORE,Number=1,Type=Float,Description="Ensembl API Scores of affinity loss/gain for TFBS measured as affinity score of alt sequence relative to the most optimal binding">';
my $info_line130 =  '##INFO=<ID=ENSTFBSANNO_TFBS7_INFO_CONTENT,Number=1,Type=String,Description="Ensembl API information content across the TFBS">';
my $info_line131 =  '##INFO=<ID=ENSTFBSANNO_TFBS7_CONS_SCORE,Number=1,Type=Float,Description="Ensembl API GERP conservation score for the entire TFBS">';
my $info_line132 =  '##INFO=<ID=ENSTFBSANNO_TFBS8_NAME,Number=1,Type=String,Description="Ensembl API Name of the overlapping TFBS">';
my $info_line133 =  '##INFO=<ID=ENSTFBSANNO_TFBS8_MATRIX,Number=1,Type=String,Description="Ensembl API Name of position weight matrix belonging to the TFBS">';
my $info_line134 =  '##INFO=<ID=ENSTFBSANNO_TFBS8_POS,Number=1,Type=Integer,Description="Ensembl API Position of variant within the TFBS">';
my $info_line135 =  '##INFO=<ID=ENSTFBSANNO_TFBS8_LENGTH,Number=1,Type=Integer,Description="Ensembl API Length of the TFBS">';
my $info_line136 =  '##INFO=<ID=ENSTFBSANNO_TFBS8_DIST_TO_START,Number=1,Type=Integer,Description="Ensembl API Distance from variant to start of TFBS">';
my $info_line137 =  '##INFO=<ID=ENSTFBSANNO_TFBS8_SIZE_OF_INDEL,Number=1,Type=Integer,Description="Ensembl API Size of the insertion or size of the deleted part of the TFBS">';
my $info_line138 =  '##INFO=<ID=ENSTFBSANNO_TFBS8_REF_SCORE,Number=1,Type=Float,Description="Ensembl API Score of affinity loss/gain for TFBS measured as affinity score of ref sequence relative to the most optimal binding">';
my $info_line139 =  '##INFO=<ID=ENSTFBSANNO_TFBS8_ALT_SCORE,Number=1,Type=Float,Description="Ensembl API Scores of affinity loss/gain for TFBS measured as affinity score of alt sequence relative to the most optimal binding">';
my $info_line140 =  '##INFO=<ID=ENSTFBSANNO_TFBS8_INFO_CONTENT,Number=1,Type=String,Description="Ensembl API information content across the TFBS">';
my $info_line141 =  '##INFO=<ID=ENSTFBSANNO_TFBS8_CONS_SCORE,Number=1,Type=Float,Description="Ensembl API GERP conservation score for the entire TFBS">';
my $info_line142 =  '##INFO=<ID=ENSTFBSANNO_TFBS9_NAME,Number=1,Type=String,Description="Ensembl API Name of the overlapping TFBS">';
my $info_line143 =  '##INFO=<ID=ENSTFBSANNO_TFBS9_MATRIX,Number=1,Type=String,Description="Ensembl API Name of position weight matrix belonging to the TFBS">';
my $info_line144 =  '##INFO=<ID=ENSTFBSANNO_TFBS9_POS,Number=1,Type=Integer,Description="Ensembl API Position of variant within the TFBS">';
my $info_line145 =  '##INFO=<ID=ENSTFBSANNO_TFBS9_LENGTH,Number=1,Type=Integer,Description="Ensembl API Length of the TFBS">';
my $info_line146 =  '##INFO=<ID=ENSTFBSANNO_TFBS9_DIST_TO_START,Number=1,Type=Integer,Description="Ensembl API Distance from variant to start of TFBS">';
my $info_line147 =  '##INFO=<ID=ENSTFBSANNO_TFBS9_SIZE_OF_INDEL,Number=1,Type=Integer,Description="Ensembl API Size of the insertion or size of the deleted part of the TFBS">';
my $info_line148 =  '##INFO=<ID=ENSTFBSANNO_TFBS9_REF_SCORE,Number=1,Type=Float,Description="Ensembl API Score of affinity loss/gain for TFBS measured as affinity score of ref sequence relative to the most optimal binding">';
my $info_line149 =  '##INFO=<ID=ENSTFBSANNO_TFBS9_ALT_SCORE,Number=1,Type=Float,Description="Ensembl API Scores of affinity loss/gain for TFBS measured as affinity score of alt sequence relative to the most optimal binding">';
my $info_line150 =  '##INFO=<ID=ENSTFBSANNO_TFBS9_INFO_CONTENT,Number=1,Type=String,Description="Ensembl API information content across the TFBS">';
my $info_line151 =  '##INFO=<ID=ENSTFBSANNO_TFBS9_CONS_SCORE,Number=1,Type=Float,Description="Ensembl API GERP conservation score for the entire TFBS">';
my $info_line152 =  '##INFO=<ID=ENSTFBSANNO_TFBS10_NAME,Number=1,Type=String,Description="Ensembl API Name of the overlapping TFBS">';
my $info_line153 =  '##INFO=<ID=ENSTFBSANNO_TFBS10_MATRIX,Number=1,Type=String,Description="Ensembl API Name of position weight matrix belonging to the TFBS">';
my $info_line154 =  '##INFO=<ID=ENSTFBSANNO_TFBS10_POS,Number=1,Type=Integer,Description="Ensembl API Position of variant within the TFBS">';
my $info_line155 =  '##INFO=<ID=ENSTFBSANNO_TFBS10_LENGTH,Number=1,Type=Integer,Description="Ensembl API Length of the TFBS">';
my $info_line156 =  '##INFO=<ID=ENSTFBSANNO_TFBS10_DIST_TO_START,Number=1,Type=Integer,Description="Ensembl API Distance from variant to start of TFBS">';
my $info_line157 =  '##INFO=<ID=ENSTFBSANNO_TFBS10_SIZE_OF_INDEL,Number=1,Type=Integer,Description="Ensembl API Size of the insertion or size of the deleted part of the TFBS">';
my $info_line158 =  '##INFO=<ID=ENSTFBSANNO_TFBS10_REF_SCORE,Number=1,Type=Float,Description="Ensembl API Score of affinity loss/gain for TFBS measured as affinity score of ref sequence relative to the most optimal binding">';
my $info_line159 =  '##INFO=<ID=ENSTFBSANNO_TFBS10_ALT_SCORE,Number=1,Type=Float,Description="Ensembl API Scores of affinity loss/gain for TFBS measured as affinity score of alt sequence relative to the most optimal binding">';
my $info_line160 =  '##INFO=<ID=ENSTFBSANNO_TFBS10_INFO_CONTENT,Number=1,Type=String,Description="Ensembl API information content across the TFBS">';
my $info_line161 =  '##INFO=<ID=ENSTFBSANNO_TFBS10_CONS_SCORE,Number=1,Type=Float,Description="Ensembl API GERP conservation score for the entire TFBS">';

# For PolyPhen and SIFT predictions
my $info_line162 =  '##INFO=<ID=ENSPREDICTANNO_POLYPHEN_TRANSCRIPT_NAME,Number=1,Type=String,Description="Ensembl API Name of transcript with the most severe PolyPhen score">';
my $info_line163 =  '##INFO=<ID=ENSPREDICTANNO_POLYPHEN_PREDICTION,Number=1,Type=String,Description="Ensembl API PolyPhen prediction for transcript">';
my $info_line164 =  '##INFO=<ID=ENSPREDICTANNO_POLYPHEN_SCORE,Number=1,Type=Float,Description="Ensembl API PolyPhen score for the transcript">';
my $info_line165 =  '##INFO=<ID=ENSPREDICTANNO_SIFT_TRANSCRIPT_NAME,Number=1,Type=String,Description="Ensembl API Name of transcript with the most severe SIFT score">';
my $info_line166 =  '##INFO=<ID=ENSPREDICTANNO_SIFT_PREDICTION,Number=1,Type=String,Description="Ensembl API SIFT prediction for transcript">';
my $info_line167 =  '##INFO=<ID=ENSPREDICTANNO_SIFT_SCORE,Number=1,Type=Float,Description="Ensembl API SIFT score for the transcript">';

###################
##### Options #####
###################

# Creating options system allowing one to decide what needs to be annotated
my ($pcanno, $lncanno, $sncanno, $reganno, $tfbsanno, $allanno, $geneint, $regint, $tfbsint, $predictanno, $isfix, $inputFile, $outputFile) = (0,0,0,0,0,0,0,0,0,0,0,0);

# Setting usage message
my $usage="\nUsage: 
	PatternAnnotator.pl -in vcf_file_input.vcf -out vcf_file_output.vcf
optional annotation settings:
	-all for all annotations
	-pc for protein coding only
	-lnc for long non coding only
	-snc for short non coding only
	-reg for regulatory only
	-tfbs for transcription factor binding sites only
	-predict for PolyPhen and SIFT predictions
	-fix if the script should run in fix mode (see documentation)
optional parameters:
	-gene_interval integer for chosing a custom interval for finding nearest gene
	-reg_interval integer for chosing a custom interval for finding nearest ENSEMBL RegulatoryFeature
	-tfbs_interval integer for chosing a custom interval for finding nearest Transcription Factor Binding Motif
	
	\nPlease try again.\n\n\n";

# Input and output are the only mandatory options
die $usage unless GetOptions(
    'pc'	=> \$pcanno,
    'lnc'	=> \$lncanno,
    'snc'	=> \$sncanno,
    'reg'	=> \$reganno,
    'tfbs'	=> \$tfbsanno,
    'predict'	=> \$predictanno,
    'all'	=> \$allanno,
    'gene_interval:i' 	=> \$geneint,
    'reg_interval:i'	=> \$regint,
    'tfbs_interval:i'	=> \$tfbsint,
    'fix' => \$isfix,
    'in:s' => \$inputFile,
    'out:s' => \$outputFile
	)
	&& $inputFile && $outputFile;

$geneint ||= 5000;
$regint ||= 5000;
$tfbsint ||= 5000;

# Opening file and creating output file
open (FILE, $inputFile) or die $!;
open (OUTPUT, "> $outputFile");

# Printing out what coordinate version is used:
my $cs = $cs_adaptor->fetch_by_name('chromosome');
print STDOUT "Coordinate System:", "\n";
printf STDOUT "\t%s %s\n", $cs->name(), $cs->version();

# Printing our what API version is used:
my $soft_version = software_version(); 
print STDOUT "ENSEMBL API Version:", "\n";
print STDOUT "\t".$soft_version."\n";

# Printing what option is used
print STDOUT "Annotations:", "\n";
if ($pcanno == 1) {
	print STDOUT "\t", "Protein Coding", "\n";
}
if ($lncanno == 1) {
	print STDOUT "\t", "Long Non-Coding", "\n";
}
if ($sncanno == 1) {
	print STDOUT "\t", "Short Non-Coding", "\n";
}
if ($reganno == 1) {
	print STDOUT "\t", "Regulatory Regions", "\n";
}
if ($tfbsanno == 1) {
	print STDOUT "\t", "Transcription Factor Binding Sites", "\n";
}
if ($predictanno == 1) {
	print STDOUT "\t", "PolyPhen and SIFT Predictions", "\n";
}
if ($allanno == 1) {
	print STDOUT "\t", "All", "\n";
}
if ($isfix == 1) {
	print STDOUT "\t", "Fix Mode", "\n";
}
if ($allanno == 1 || $pcanno == 1 || $lncanno == 1 || $sncanno == 1 || $reganno == 1) {
	print STDOUT "Intervals Around Variants for Finding Closest Features:", "\n";
}
if ($allanno == 1 || $pcanno == 1 || $lncanno == 1 || $sncanno == 1) {
	print STDOUT "\t"."Genes: ".$geneint." bases"."\n";
}
if ($allanno == 1 || $reganno == 1 || $tfbsanno == 1) {
	print STDOUT "\t"."Regulatory Features: ".$regint." bases"."\n";
	print STDOUT "\t"."TFBS: ".$tfbsint." bases"."\n";
}
print STDOUT "\n";

###############################################
##### Going through the file line by line #####
###############################################

while (<FILE>) {
    chomp;
    my $line = $_;

	# Adding appropriate header information
    if ($isfix == 0) {
		if ($line =~ /^#/) {
			print OUTPUT $line, "\n";
			if ($pcanno == 1 || $allanno == 1) {  
				if ($line =~ /^##INFO=<ID=END/) {
					print OUTPUT $info_line1, "\n", $info_line2, "\n", $info_line3, "\n", $info_line4, "\n", $info_line5, "\n", $info_line6, "\n", $info_line7, "\n", $info_line8, "\n", $info_line9, "\n", $info_line10, "\n", $info_line11, "\n", $info_line12, "\n";
				}
			}
			if ($lncanno == 1 || $allanno == 1) {  
				if ($line =~ /^##INFO=<ID=END/) {
					print OUTPUT $info_line13, "\n", $info_line14, "\n", $info_line15, "\n", $info_line16, "\n", $info_line17, "\n", $info_line18, "\n", $info_line19, "\n", $info_line20, "\n", $info_line21, "\n", $info_line22, "\n", $info_line23, "\n";
				}
			}
			if ($sncanno == 1 || $allanno == 1) {  
				if ($line =~ /^##INFO=<ID=END/) {
					print OUTPUT $info_line24, "\n", $info_line25, "\n", $info_line26, "\n", $info_line27, "\n", $info_line28, "\n", $info_line29, "\n", $info_line30, "\n", $info_line31, "\n", $info_line32, "\n", $info_line33, "\n", $info_line34, "\n";
				}
			}
			if ($reganno == 1 || $allanno == 1) {  
				if ($line =~ /^##INFO=<ID=END/) {
					print OUTPUT $info_line35, "\n", $info_line36, "\n", $info_line37, "\n", $info_line38, "\n", $info_line39, "\n", $info_line40, "\n", $info_line41, "\n", $info_line42, "\n", $info_line43, "\n", $info_line44, "\n", $info_line45, "\n", $info_line46, "\n", $info_line47, "\n", $info_line48, "\n", $info_line49, "\n", $info_line50, "\n", $info_line51, "\n", $info_line52, "\n", $info_line53, "\n", $info_line54, "\n", $info_line55, "\n", $info_line56, "\n", $info_line57, "\n", $info_line58, "\n", $info_line59, "\n", $info_line60, "\n", $info_line61, "\n";
				}
			}
			if ($tfbsanno == 1 || $allanno == 1) {  
				if ($line =~ /^##INFO=<ID=END/) {
					print OUTPUT $info_line62, "\n", $info_line63, "\n", $info_line64, "\n", $info_line65, "\n", $info_line66, "\n", $info_line67, "\n", $info_line68, "\n", $info_line69, "\n", $info_line70, "\n", $info_line71, "\n", $info_line72, "\n", $info_line73, "\n", $info_line74, "\n", $info_line75, "\n", $info_line76, "\n", $info_line77, "\n", $info_line78, "\n", $info_line79, "\n", $info_line80, "\n", $info_line81, "\n", $info_line82, "\n", $info_line83, "\n", $info_line84, "\n", $info_line85, "\n", $info_line86, "\n", $info_line87, "\n", $info_line88, "\n", $info_line89, "\n", $info_line90, "\n", $info_line91, "\n", $info_line92, "\n", $info_line93, "\n", $info_line94, "\n", $info_line95, "\n", $info_line96, "\n", $info_line97, "\n", $info_line98, "\n", $info_line99, "\n", $info_line100, "\n", $info_line101, "\n", $info_line102, "\n", $info_line103, "\n", $info_line104, "\n", $info_line105, "\n", $info_line106, "\n", $info_line107, "\n", $info_line108, "\n", $info_line109, "\n", $info_line110, "\n", $info_line111, "\n", $info_line112, "\n", $info_line113, "\n", $info_line114, "\n", $info_line115, "\n", $info_line116, "\n", $info_line117, "\n", $info_line118, "\n", $info_line119, "\n", $info_line120, "\n", $info_line121, "\n", $info_line122, "\n", $info_line123, "\n", $info_line124, "\n", $info_line125, "\n", $info_line126, "\n", $info_line127, "\n", $info_line128, "\n", $info_line129, "\n", $info_line130, "\n", $info_line131, "\n", $info_line132, "\n", $info_line133, "\n", $info_line134, "\n", $info_line135, "\n", $info_line136, "\n", $info_line137, "\n", $info_line138, "\n", $info_line139, "\n", $info_line140, "\n", $info_line141, "\n", $info_line142, "\n", $info_line143, "\n", $info_line144, "\n", $info_line145, "\n", $info_line146, "\n", $info_line147, "\n", $info_line148, "\n", $info_line149, "\n", $info_line150, "\n", $info_line151, "\n", $info_line152, "\n", $info_line153, "\n", $info_line154, "\n", $info_line155, "\n", $info_line156, "\n", $info_line157, "\n", $info_line158, "\n", $info_line159, "\n", $info_line160, "\n", $info_line161, "\n";
				}
			}
			if ($predictanno == 1 || $allanno == 1) {  
				if ($line =~ /^##INFO=<ID=END/) {
					print OUTPUT $info_line162, "\n", $info_line163, "\n", $info_line164, "\n", $info_line165, "\n", $info_line166, "\n", $info_line167, "\n";
				}
			}
			
			next;
		}
	}
	
	# In fix mode the header from the original vcf-file is printed as it is
	if ($isfix == 1) {
		if ($line =~ /^#/) {
			print OUTPUT $line, "\n";
		}
	}
    
    ############################################################
    ##### Cleaning variant and creating a VariationFeature #####
    ############################################################
    
    # Container for annotation string
    my $annotationString ="";
    
    my $vf = ();
    if ($line !~ /^#/) {
		my @col = split("\t", $line);
		my $slice = $slice_adaptor->fetch_by_region('chromosome', $col[0]);
		
		# A VariationFeature is created for variants where at least one of the alleles is reference
    	if ($col[4] !~ /,/) {
    		# Clean SNPs
    		my @cleaningInfoSNP = PatternModule::CleaningVariants::cleanSNP($col[3], $col[4]);
    		if (scalar @cleaningInfoSNP == 3) {
    			print STDOUT "Cleaning complex variant:", "\n";
   	 			print STDOUT "\tFrom: " . $col[3] . "/" . $col[4] . "\n";
   	 			print STDOUT "\tStart position: " . $col[1] . "\n";
   	 			$col[3] = $cleaningInfoSNP[0];
   	 			$col[4] = $cleaningInfoSNP[1];
   	 			$col[1] = $col[1] + $cleaningInfoSNP[2];
   	 			print STDOUT "\tTo: " . $col[3] . "/" . $col[4] . "\n";
   	 			print STDOUT "\t" . "Start position: " . $col[1] . "\n";
   	 		}
   	 		# Clean Deletions
   	 		my @cleaningInfoDELETION = PatternModule::CleaningVariants::cleanDELETIONS($col[3], $col[4]);
   	 		if (scalar @cleaningInfoDELETION == 3) {
    			print STDOUT "Cleaning deletion variant:", "\n";
   	 			print STDOUT "\tFrom: " . $col[3] . "/" . $col[4] . "\n";
   	 			print STDOUT "\tStart position: " . $col[1] . "\n";
   	 			$col[3] = $cleaningInfoDELETION[0];
   	 			$col[4] = $cleaningInfoDELETION[1]; 
   	 			$col[1] = $col[1] + $cleaningInfoDELETION[2];	
   	 			print STDOUT "\tTo: " . $col[3] . "/" . $col[4] . "\n";	
   	 			print STDOUT "\t" . "Start position: " . $col[1] . "\n";
   	 		}
   	 		# Clean Insertions
   	 		my @cleaningInfoINSERTION = PatternModule::CleaningVariants::cleanINSERTIONS($col[3], $col[4]);
   	 		if (scalar @cleaningInfoINSERTION == 3) {
    			print STDOUT "Cleaning insertion variant:", "\n";
   	 			print STDOUT "\tFrom: " . $col[3] . "/" . $col[4] . "\n";
   	 			print STDOUT "\tStart position: " . $col[1] . "\n";
   	 			$col[3] = $cleaningInfoINSERTION[0];
   	 			$col[4] = $cleaningInfoINSERTION[1]; 
   	 			$col[1] = $col[1] + $cleaningInfoINSERTION[2];	
   	 			print STDOUT "\tTo: " . $col[3] . "/" . $col[4] . "\n";	
   	 			print STDOUT "\t" . "Start position: " . $col[1] . "\n";
   	 		}
    		$vf = PatternModule::GeneAnno::createVariationFeature($col[1], $col[3], $col[4], $slice, $variant_feature_adaptor);
    	}
    	
    	# Creating one VariationFeature for each multi allelic variant
    	if ($col[4] =~ /,/) {
    		$vf = PatternModule::GeneAnno::createMultiAllelicVariationFeature($col[1], $col[3], $col[4], $slice, $variant_feature_adaptor);
    	}
    	
		####################
		##### Fix Mode #####
		####################
	
		# Lines should be annotated per default unless in fix mode
		my $annotate = 1;
		
		# In fix mode no lines are re-annotated per default
		if ($isfix == 1) {
			$annotate = 0;
		}
		
		# Fix mode is activated only for lines containing the information IS_MULTI_ALLELIC=YES
		if ($line =~ /IS_MULTI_ALLELIC=YES/) {
			$annotate = 1;
		}
    	
    	################################################
    	##### Annotation of protein coding regions #####
    	################################################
    	
    	if ($annotate == 1) {
			if ($allanno == 1 || $pcanno == 1) {
				my %variantInformation = PatternModule::GeneAnno::createVariantProteinCodingHash();
				if ($vf) {
					my ($one, $two) = PatternModule::GeneAnno::doesVariantOverlapWithProteinCodingGene($vf);
					$variantInformation{'ENSPCANNO_VARIANT_OVERLAPS_WITH_GENE'} = $one;
					$variantInformation{'ENSPCANNO_GENE'} = $two;
					my ($three, $four) = PatternModule::GeneAnno::findMostSeverelyHitTranscript($vf,$transcript_adaptor, "PC");
					$variantInformation{'ENSPCANNO_TRANSCRIPT'} = $three;
					$variantInformation{'ENSPCANNO_TRANSCRIPT_CONSEQUENCE'} = $four;
					if ($three ne "NA") {
						$variantInformation{'ENSPCANNO_GENE'} = $transcript_adaptor->fetch_by_stable_id($three)->get_Gene()->display_xref->display_id();
						$variantInformation{'ENSPCANNO_TRANSCRIPT_LENGTH'} = $transcript_adaptor->fetch_by_stable_id($three)->seq_region_end() - $transcript_adaptor->fetch_by_stable_id($three)->seq_region_start() + 1;
					}
					if ($variantInformation{'ENSPCANNO_GENE'} eq "NA") {
						my ($five, $six) = PatternModule::GeneAnno::getNearestGene($vf, $slice_adaptor, "PC", $geneint);
						$variantInformation{'ENSPCANNO_GENE'} = $five;
						$variantInformation{'ENSPCANNO_TRANSCRIPT_CONSEQUENCE'} = $six;
					}
					if ($variantInformation{'ENSPCANNO_GENE'} ne "NA") {
						my $gene = $gene_adaptor->fetch_by_display_label($variantInformation{'ENSPCANNO_GENE'});
						$variantInformation{'ENSPCANNO_GENE_LENGTH'} = $gene->length();	
						$variantInformation{'ENSPCANNO_GENE_NUMBER_OF_TRANSCRIPTS'} = scalar @{$gene->get_all_Transcripts()};				 
						my ($seven, $eight) = PatternModule::GeneAnno::findDistanceToStartAndEndOfGene($gene, $vf);
						$variantInformation{'ENSPCANNO_GENE_DISTANCE_TO_START'} = $seven;
						$variantInformation{'ENSPCANNO_GENE_DISTANCE_TO_END'} = $eight;
					}
					if ($variantInformation{'ENSPCANNO_TRANSCRIPT'} ne "NA") {
						my $transcript = $transcript_adaptor->fetch_by_stable_id($variantInformation{'ENSPCANNO_TRANSCRIPT'});
						my ($nine, $ten) = PatternModule::GeneAnno::findDistanceToStartAndEndOfTranscript($transcript, $vf);
						$variantInformation{'ENSPCANNO_TRANSCRIPT_DISTANCE_TO_START'} = $nine;
						$variantInformation{'ENSPCANNO_TRANSCRIPT_DISTANCE_TO_END'} = $ten;
						$variantInformation{'ENSPCANNO_TRANSCRIPT_BIOTYPE'} = "protein_coding";
					}
				}
			
				# Printing annotation information to terminal
				print STDOUT $variantInformation{'ENSPCANNO_VARIANT_OVERLAPS_WITH_GENE'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSPCANNO_GENE'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSPCANNO_GENE_LENGTH'};
				print STDOUT "\t";				
				print STDOUT $variantInformation{'ENSPCANNO_GENE_NUMBER_OF_TRANSCRIPTS'};
				print STDOUT "\t";		
				print STDOUT $variantInformation{'ENSPCANNO_GENE_DISTANCE_TO_START'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSPCANNO_GENE_DISTANCE_TO_END'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSPCANNO_TRANSCRIPT'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSPCANNO_TRANSCRIPT_BIOTYPE'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSPCANNO_TRANSCRIPT_LENGTH'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSPCANNO_TRANSCRIPT_DISTANCE_TO_START'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSPCANNO_TRANSCRIPT_DISTANCE_TO_END'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSPCANNO_TRANSCRIPT_CONSEQUENCE'};
				print STDOUT "\t";
			
				# Adding information to the output vcf-file
				$annotationString .= 
					";ENSPCANNO_VARIANT_OVERLAPS_WITH_GENE=" . $variantInformation{'ENSPCANNO_VARIANT_OVERLAPS_WITH_GENE'} .
					";ENSPCANNO_GENE=" . $variantInformation{'ENSPCANNO_GENE'} .
					";ENSPCANNO_GENE_LENGTH=" . $variantInformation{'ENSPCANNO_GENE_LENGTH'} .
					";ENSPCANNO_GENE_NUMBER_OF_TRANSCRIPTS=" . $variantInformation{'ENSPCANNO_GENE_NUMBER_OF_TRANSCRIPTS'} .
					";ENSPCANNO_GENE_DISTANCE_TO_START=" . $variantInformation{'ENSPCANNO_GENE_DISTANCE_TO_START'} .
					";ENSPCANNO_GENE_DISTANCE_TO_END=" . $variantInformation{'ENSPCANNO_GENE_DISTANCE_TO_END'} .
					";ENSPCANNO_TRANSCRIPT=" . $variantInformation{'ENSPCANNO_TRANSCRIPT'} .
					";ENSPCANNO_TRANSCRIPT_BIOTYPE=" . $variantInformation{'ENSPCANNO_TRANSCRIPT_BIOTYPE'} .
					";ENSPCANNO_TRANSCRIPT_LENGTH=" . $variantInformation{'ENSPCANNO_TRANSCRIPT_LENGTH'} .
					";ENSPCANNO_TRANSCRIPT_DISTANCE_TO_START=" . $variantInformation{'ENSPCANNO_TRANSCRIPT_DISTANCE_TO_START'} .
					";ENSPCANNO_TRANSCRIPT_DISTANCE_TO_END=" . $variantInformation{'ENSPCANNO_TRANSCRIPT_DISTANCE_TO_END'} .
					";ENSPCANNO_TRANSCRIPT_CONSEQUENCE=" . $variantInformation{'ENSPCANNO_TRANSCRIPT_CONSEQUENCE'};	
			}
		
			#################################################
			##### Annotation of long non-coding regions #####
			#################################################
		
			if ($allanno == 1 || $lncanno == 1) {
				my %variantInformation = PatternModule::GeneAnno::createVariantLongNonCodingHash();
				if ($vf) {
					my ($one) = PatternModule::GeneAnno::doesVariantOverlapWithLongNonCodingGene($vf);
					$variantInformation{'ENSLNCANNO_GENE'} = $one;				
					if ($variantInformation{'ENSLNCANNO_GENE'} eq "NA") {
						my ($two, $three) = PatternModule::GeneAnno::getNearestGene($vf, $slice_adaptor, "LNC", $geneint);
						$variantInformation{'ENSLNCANNO_GENE'} = $two;
						$variantInformation{'ENSLNCANNO_TRANSCRIPT_CONSEQUENCE'} = $three;
					}
					my ($four, $five) = PatternModule::GeneAnno::findMostSeverelyHitTranscript($vf,$transcript_adaptor, "LNC");
					$variantInformation{'ENSLNCANNO_TRANSCRIPT'} = $four;
					if ($five ne "NA") {
						$variantInformation{'ENSLNCANNO_TRANSCRIPT_CONSEQUENCE'} = $five;
					}
					if ($four ne "NA") {
						$variantInformation{'ENSLNCANNO_GENE'} = $transcript_adaptor->fetch_by_stable_id($four)->get_Gene()->display_xref->display_id();
						$variantInformation{'ENSLNCANNO_TRANSCRIPT_LENGTH'} = $transcript_adaptor->fetch_by_stable_id($four)->seq_region_end() - $transcript_adaptor->fetch_by_stable_id($four)->seq_region_start() + 1;
					}
					if ($variantInformation{'ENSLNCANNO_GENE'} ne "NA") {
						my $gene = $gene_adaptor->fetch_by_display_label($variantInformation{'ENSLNCANNO_GENE'});
						$variantInformation{'ENSLNCANNO_GENE_LENGTH'} = $gene->length();	
						$variantInformation{'ENSLNCANNO_GENE_NUMBER_OF_TRANSCRIPTS'} = scalar @{$gene->get_all_Transcripts()};				 
						my ($six, $seven) = PatternModule::GeneAnno::findDistanceToStartAndEndOfGene($gene, $vf);
						$variantInformation{'ENSLNCANNO_GENE_DISTANCE_TO_START'} = $six;
						$variantInformation{'ENSLNCANNO_GENE_DISTANCE_TO_END'} = $seven;
						if ($variantInformation{'ENSLNCANNO_TRANSCRIPT'} ne "NA") {
							my $transcript = $transcript_adaptor->fetch_by_stable_id($variantInformation{'ENSLNCANNO_TRANSCRIPT'});
							$variantInformation{'ENSLNCANNO_TRANSCRIPT_BIOTYPE'} = $transcript->biotype;
							my ($eight, $nine) = PatternModule::GeneAnno::findDistanceToStartAndEndOfTranscript($transcript, $vf);
							$variantInformation{'ENSLNCANNO_TRANSCRIPT_DISTANCE_TO_START'} = $eight;
							$variantInformation{'ENSLNCANNO_TRANSCRIPT_DISTANCE_TO_END'} = $nine;
						}
					}
				
				}
			
				# Printing annotation information to terminal
				print STDOUT $variantInformation{'ENSLNCANNO_GENE'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSLNCANNO_GENE_LENGTH'};
				print STDOUT "\t";				
				print STDOUT $variantInformation{'ENSLNCANNO_GENE_NUMBER_OF_TRANSCRIPTS'};
				print STDOUT "\t";		
				print STDOUT $variantInformation{'ENSLNCANNO_GENE_DISTANCE_TO_START'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSLNCANNO_GENE_DISTANCE_TO_END'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSLNCANNO_TRANSCRIPT'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSLNCANNO_TRANSCRIPT_BIOTYPE'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSLNCANNO_TRANSCRIPT_LENGTH'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSLNCANNO_TRANSCRIPT_DISTANCE_TO_START'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSLNCANNO_TRANSCRIPT_DISTANCE_TO_END'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSLNCANNO_TRANSCRIPT_CONSEQUENCE'};
				print STDOUT "\t";

			
				# Adding information to the output vcf-file
				$annotationString .= 
					";ENSLNCANNO_GENE=" . $variantInformation{'ENSLNCANNO_GENE'} .
					";ENSLNCANNO_GENE_LENGTH=" . $variantInformation{'ENSLNCANNO_GENE_LENGTH'} .
					";ENSLNCANNO_GENE_NUMBER_OF_TRANSCRIPTS=" . $variantInformation{'ENSLNCANNO_GENE_NUMBER_OF_TRANSCRIPTS'} .
					";ENSLNCANNO_GENE_DISTANCE_TO_START=" . $variantInformation{'ENSLNCANNO_GENE_DISTANCE_TO_START'} .
					";ENSLNCANNO_GENE_DISTANCE_TO_END=" . $variantInformation{'ENSLNCANNO_GENE_DISTANCE_TO_END'} .
					";ENSLNCANNO_TRANSCRIPT=" . $variantInformation{'ENSLNCANNO_TRANSCRIPT'} .
					";ENSLNCANNO_TRANSCRIPT_BIOTYPE=" . $variantInformation{'ENSLNCANNO_TRANSCRIPT_BIOTYPE'} .
					";ENSLNCANNO_TRANSCRIPT_LENGTH=" . $variantInformation{'ENSLNCANNO_TRANSCRIPT_LENGTH'} .
					";ENSLNCANNO_TRANSCRIPT_DISTANCE_TO_START=" . $variantInformation{'ENSLNCANNO_TRANSCRIPT_DISTANCE_TO_START'} .
					";ENSLNCANNO_TRANSCRIPT_DISTANCE_TO_END=" . $variantInformation{'ENSLNCANNO_TRANSCRIPT_DISTANCE_TO_END'} .
					";ENSLNCANNO_TRANSCRIPT_CONSEQUENCE=" . $variantInformation{'ENSLNCANNO_TRANSCRIPT_CONSEQUENCE'};	
			}
		
			##################################################
			##### Annotation of short non-coding regions #####
			##################################################
		
			if ($allanno == 1 || $sncanno == 1) {
				my %variantInformation = PatternModule::GeneAnno::createVariantShortNonCodingHash();
				if ($vf) {
					my ($one) = PatternModule::GeneAnno::doesVariantOverlapWithShortNonCodingGene($vf);
					$variantInformation{'ENSSNCANNO_GENE'} = $one;				
					if ($variantInformation{'ENSSNCANNO_GENE'} eq "NA") {
						my ($two, $three) = PatternModule::GeneAnno::getNearestGene($vf, $slice_adaptor, "SNC", $geneint);
						$variantInformation{'ENSSNCANNO_GENE'} = $two;
						$variantInformation{'ENSSNCANNO_TRANSCRIPT_CONSEQUENCE'} = $three;
					}
					my ($four, $five) = PatternModule::GeneAnno::findMostSeverelyHitTranscript($vf,$transcript_adaptor, "SNC");
					$variantInformation{'ENSSNCANNO_TRANSCRIPT'} = $four;
					if ($five ne "NA") {
						$variantInformation{'ENSSNCANNO_TRANSCRIPT_CONSEQUENCE'} = $five;
					}
					if ($four ne "NA") {
						$variantInformation{'ENSSNCANNO_GENE'} = $transcript_adaptor->fetch_by_stable_id($four)->get_Gene()->display_xref->display_id();
						$variantInformation{'ENSSNCANNO_TRANSCRIPT_LENGTH'} = $transcript_adaptor->fetch_by_stable_id($four)->seq_region_end() - $transcript_adaptor->fetch_by_stable_id($four)->seq_region_start() + 1;
					}
					if ($variantInformation{'ENSSNCANNO_GENE'} ne "NA") {
						my $gene = $gene_adaptor->fetch_by_display_label($variantInformation{'ENSSNCANNO_GENE'});
						$variantInformation{'ENSSNCANNO_GENE_LENGTH'} = $gene->length();	
						$variantInformation{'ENSSNCANNO_GENE_NUMBER_OF_TRANSCRIPTS'} = scalar @{$gene->get_all_Transcripts()};				 
						my ($six, $seven) = PatternModule::GeneAnno::findDistanceToStartAndEndOfGene($gene, $vf);
						$variantInformation{'ENSSNCANNO_GENE_DISTANCE_TO_START'} = $six;
						$variantInformation{'ENSSNCANNO_GENE_DISTANCE_TO_END'} = $seven;
						if ($variantInformation{'ENSSNCANNO_TRANSCRIPT'} ne "NA") {
							my $transcript = $transcript_adaptor->fetch_by_stable_id($variantInformation{'ENSSNCANNO_TRANSCRIPT'});
							$variantInformation{'ENSSNCANNO_TRANSCRIPT_BIOTYPE'} = $transcript->biotype();
							my ($eight, $nine) = PatternModule::GeneAnno::findDistanceToStartAndEndOfTranscript($transcript, $vf);
							$variantInformation{'ENSSNCANNO_TRANSCRIPT_DISTANCE_TO_START'} = $eight;
							$variantInformation{'ENSSNCANNO_TRANSCRIPT_DISTANCE_TO_END'} = $nine;
						}
					}
				}
			
				# Printing annotation information to terminal
				print STDOUT $variantInformation{'ENSSNCANNO_GENE'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSSNCANNO_GENE_LENGTH'};
				print STDOUT "\t";				
				print STDOUT $variantInformation{'ENSSNCANNO_GENE_NUMBER_OF_TRANSCRIPTS'};
				print STDOUT "\t";		
				print STDOUT $variantInformation{'ENSSNCANNO_GENE_DISTANCE_TO_START'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSSNCANNO_GENE_DISTANCE_TO_END'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSSNCANNO_TRANSCRIPT'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSSNCANNO_TRANSCRIPT_BIOTYPE'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSSNCANNO_TRANSCRIPT_LENGTH'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSSNCANNO_TRANSCRIPT_DISTANCE_TO_START'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSSNCANNO_TRANSCRIPT_DISTANCE_TO_END'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSSNCANNO_TRANSCRIPT_CONSEQUENCE'};
				print STDOUT "\t";

				# Adding information to the output vcf-file
				$annotationString .= 
					";ENSSNCANNO_GENE=" . $variantInformation{'ENSSNCANNO_GENE'} .
					";ENSSNCANNO_GENE_LENGTH=" . $variantInformation{'ENSSNCANNO_GENE_LENGTH'} .
					";ENSSNCANNO_GENE_NUMBER_OF_TRANSCRIPTS=" . $variantInformation{'ENSSNCANNO_GENE_NUMBER_OF_TRANSCRIPTS'} .
					";ENSSNCANNO_GENE_DISTANCE_TO_START=" . $variantInformation{'ENSSNCANNO_GENE_DISTANCE_TO_START'} .
					";ENSSNCANNO_GENE_DISTANCE_TO_END=" . $variantInformation{'ENSSNCANNO_GENE_DISTANCE_TO_END'} .
					";ENSSNCANNO_TRANSCRIPT=" . $variantInformation{'ENSSNCANNO_TRANSCRIPT'} .
					";ENSSNCANNO_TRANSCRIPT_BIOTYPE=" . $variantInformation{'ENSSNCANNO_TRANSCRIPT_BIOTYPE'} .
					";ENSSNCANNO_TRANSCRIPT_LENGTH=" . $variantInformation{'ENSSNCANNO_TRANSCRIPT_LENGTH'} .
					";ENSSNCANNO_TRANSCRIPT_DISTANCE_TO_START=" . $variantInformation{'ENSSNCANNO_TRANSCRIPT_DISTANCE_TO_START'} .
					";ENSSNCANNO_TRANSCRIPT_DISTANCE_TO_END=" . $variantInformation{'ENSSNCANNO_TRANSCRIPT_DISTANCE_TO_END'} .
					";ENSSNCANNO_TRANSCRIPT_CONSEQUENCE=" . $variantInformation{'ENSSNCANNO_TRANSCRIPT_CONSEQUENCE'};	
			}
		
			############################################
			##### Annotation of regulatory regions #####
			############################################
		
			if ($allanno == 1 || $reganno == 1) {
				my %variantInformation = PatternModule::RegAnno::createVariantRegulatoryHash();
			
				# Checking if the variation feature overlaps with a regulatory feature
				if ($vf) {
					my ($one, $two, $three) = PatternModule::RegAnno::doesVariantOverlapWithRegulatoryFeature($vf);
					$variantInformation{'ENSREGANNO_VARIANT_OVERLAPS_WITH_REGULATORY'} = $one;
					$variantInformation{'ENSREGANNO_REGULATORY_FEATURE'} = $two;
					$variantInformation{'ENSREGANNO_CONSEQUENCE'} = $three;
				}
				
				# If the variant feature does overlap with a regulatory feature, further annotations are added
				if ($variantInformation{'ENSREGANNO_REGULATORY_FEATURE'} ne "NA") {
					my $reg_feature = $regfeat_adaptor->fetch_by_stable_id($variantInformation{'ENSREGANNO_REGULATORY_FEATURE'});
					my ($six, $seven) = PatternModule::RegAnno::findDistanceToStartAndEndOfRegFeature($reg_feature, $vf);
					$variantInformation{'ENSREGANNO_DISTANCE_TO_START'} = $six;
					$variantInformation{'ENSREGANNO_DISTANCE_TO_END'} = $seven;
					
					# Cell specific information is pulled from the regulatory feature
					my @per_cell_reg_features = @{$regfeat_adaptor->fetch_all_by_stable_ID($reg_feature->stable_id)};
				
					# Testing whether the variationfeature is multiallelic		
					my $x = "/";
					my $y = $vf->allele_string;	
					my @c = $y =~ /$x/g;
					my $count = @c;
				
					# Variants with multiple alternative alleles
					if ($count > 1) {
						my @variationFeatureList = PatternModule::GeneAnno::splitMultiAllelicVariationFeature($vf, $slice, $variant_feature_adaptor);
					
						my %multiTissueHash = ();
						my %NeuroHash = ();
						my %ImmunoHash = ();
						
						# For each of the variation features the most severely hit transcription factor binding motif is found and is placed in a hash
						foreach my $multiVf (@variationFeatureList) {
							my @multiTissueScores = ();
							my @NeuroScores = ();
							my @ImmunoScores = ();
						
							my @multiTissueInfo = PatternModule::RegAnno::getTranscriptionFactorBindingScore(\@per_cell_reg_features, $multiVf, $slice_adaptor, "Multiple_Tissues");
							my @NeuroInfo = PatternModule::RegAnno::getTranscriptionFactorBindingScore(\@per_cell_reg_features, $multiVf, $slice_adaptor, "Neuro");
							my @ImmunoInfo = PatternModule::RegAnno::getTranscriptionFactorBindingScore(\@per_cell_reg_features, $multiVf, $slice_adaptor, "Immuno");
						
							if ($multiTissueInfo[2] ne ".") {
								push @multiTissueScores, $multiTissueInfo[2];
							}
							if ($multiTissueInfo[3] ne ".") {
								push @multiTissueScores, split(",", $multiTissueInfo[3]);
							}
							if ($NeuroInfo[2] ne ".") {
								push @NeuroScores, $NeuroInfo[2];
							}
							if ($NeuroInfo[3] ne ".") {
								push @NeuroScores, split(",", $NeuroInfo[3]);
							}
							if ($ImmunoInfo[2] ne ".") {
								push @ImmunoScores, $ImmunoInfo[2];
							}
							if ($ImmunoInfo[3] ne ".") {
								push @ImmunoScores, split(",", $ImmunoInfo[3]);
							}
							if (scalar @multiTissueScores > 0) {
								$multiTissueHash{min(@multiTissueScores)} = $multiVf;
							}
							if (scalar @NeuroScores > 0) {
								$NeuroHash{min(@NeuroScores)} = $multiVf;
							}
							if (scalar @ImmunoScores > 0) {
								$ImmunoHash{min(@ImmunoScores)} = $multiVf;
							}
						}
						
						# For variation feature with the most severely hit transcription factor binding motif is used for annotating the effect of the multiallelic variant
						if(%multiTissueHash) { 
							my $multiTissueMostSevere = (sort {$a <=> $b} keys %multiTissueHash)[0];
							my ($eight, $nine, $ten, $eleven) = PatternModule::RegAnno::getTranscriptionFactorBindingScore(\@per_cell_reg_features, $multiTissueHash{$multiTissueMostSevere}, $slice_adaptor, "Multiple_Tissues");
							$variantInformation{'ENSREGANNO_MULTICELL_TFBS'} = $eight;
							$variantInformation{'ENSREGANNO_MULTICELL_MULTITFBS'} = $nine;						
							$variantInformation{'ENSREGANNO_MULTICELL_TFBS_SCORE'} = $ten;
							$variantInformation{'ENSREGANNO_MULTICELL_MULTITFBS_SCORE'} = $eleven;
						}
						else {
							my ($eight, $nine, $ten, $eleven) = PatternModule::RegAnno::getTranscriptionFactorBindingScore(\@per_cell_reg_features, $variationFeatureList[0], $slice_adaptor, "Multiple_Tissues");
							$variantInformation{'ENSREGANNO_MULTICELL_TFBS'} = $eight;
							$variantInformation{'ENSREGANNO_MULTICELL_MULTITFBS'} = $nine;
							$variantInformation{'ENSREGANNO_MULTICELL_TFBS_SCORE'} = $ten;
							$variantInformation{'ENSREGANNO_MULTICELL_MULTITFBS_SCORE'} = $eleven;
						}
						if(%NeuroHash) { 
							my $NeuroMostSevere = (sort {$a <=> $b} keys %NeuroHash)[0];
							my ($twelve, $thirteen, $fourteen, $fifteen) = PatternModule::RegAnno::getTranscriptionFactorBindingScore(\@per_cell_reg_features, $NeuroHash{$NeuroMostSevere}, $slice_adaptor, "Neuro");
							$variantInformation{'ENSREGANNO_NEUROCELL_TFBS'} = $twelve;
							$variantInformation{'ENSREGANNO_NEUROCELL_MULTITFBS'} = $thirteen;
							$variantInformation{'ENSREGANNO_NEUROCELL_TFBS_SCORE'} = $fourteen;
							$variantInformation{'ENSREGANNO_NEUROCELL_MULTITFBS_SCORE'} = $fifteen;
						}
						else {
							my ($twelve, $thirteen, $fourteen, $fifteen) = PatternModule::RegAnno::getTranscriptionFactorBindingScore(\@per_cell_reg_features, $variationFeatureList[0], $slice_adaptor, "Neuro");
							$variantInformation{'ENSREGANNO_NEUROCELL_TFBS'} = $twelve;
							$variantInformation{'ENSREGANNO_NEUROCELL_MULTITFBS'} = $thirteen;
							$variantInformation{'ENSREGANNO_NEUROCELL_TFBS_SCORE'} = $fourteen;
							$variantInformation{'ENSREGANNO_NEUROCELL_MULTITFBS_SCORE'} = $fifteen;
						}
						if(%ImmunoHash) { 
							my $ImmunoMostSevere = (sort {$a <=> $b} keys %ImmunoHash)[0];
							my ($sixteen, $seventeen, $eightteen, $nineteen) = PatternModule::RegAnno::getTranscriptionFactorBindingScore(\@per_cell_reg_features, $ImmunoHash{$ImmunoMostSevere}, $slice_adaptor, "Immuno");
							$variantInformation{'ENSREGANNO_IMMUNOCELL_TFBS'} = $sixteen;
							$variantInformation{'ENSREGANNO_IMMUNOCELL_MULTITFBS'} = $seventeen;
							$variantInformation{'ENSREGANNO_IMMUNOCELL_TFBS_SCORE'} = $eightteen;
							$variantInformation{'ENSREGANNO_IMMUNOCELL_MULTITFBS_SCORE'} = $nineteen;
						}
						else {
							my ($sixteen, $seventeen, $eightteen, $nineteen) = PatternModule::RegAnno::getTranscriptionFactorBindingScore(\@per_cell_reg_features, $variationFeatureList[0], $slice_adaptor, "Immuno");
							$variantInformation{'ENSREGANNO_IMMUNOCELL_TFBS'} = $sixteen;
							$variantInformation{'ENSREGANNO_IMMUNOCELL_MULTITFBS'} = $seventeen;
							$variantInformation{'ENSREGANNO_IMMUNOCELL_TFBS_SCORE'} = $eightteen;
							$variantInformation{'ENSREGANNO_IMMUNOCELL_MULTITFBS_SCORE'} = $nineteen;
						}	
					}
				
					# Variants with a single alternative allele
					if ($count == 1) {
						my ($eight, $nine, $ten, $eleven) = PatternModule::RegAnno::getTranscriptionFactorBindingScore(\@per_cell_reg_features, $vf, $slice_adaptor, "Multiple_Tissues");
						my ($twelve, $thirteen, $fourteen, $fifteen) = PatternModule::RegAnno::getTranscriptionFactorBindingScore(\@per_cell_reg_features, $vf, $slice_adaptor, "Neuro");
						my ($sixteen, $seventeen, $eightteen, $nineteen) = PatternModule::RegAnno::getTranscriptionFactorBindingScore(\@per_cell_reg_features, $vf, $slice_adaptor, "Immuno");
						$variantInformation{'ENSREGANNO_MULTICELL_TFBS'} = $eight;
						$variantInformation{'ENSREGANNO_MULTICELL_MULTITFBS'} = $nine;
						$variantInformation{'ENSREGANNO_MULTICELL_TFBS_SCORE'} = $ten;
						$variantInformation{'ENSREGANNO_MULTICELL_MULTITFBS_SCORE'} = $eleven;
						$variantInformation{'ENSREGANNO_NEUROCELL_TFBS'} = $twelve;
						$variantInformation{'ENSREGANNO_NEUROCELL_MULTITFBS'} = $thirteen;
						$variantInformation{'ENSREGANNO_NEUROCELL_TFBS_SCORE'} = $fourteen;
						$variantInformation{'ENSREGANNO_NEUROCELL_MULTITFBS_SCORE'} = $fifteen;
						$variantInformation{'ENSREGANNO_IMMUNOCELL_TFBS'} = $sixteen;
						$variantInformation{'ENSREGANNO_IMMUNOCELL_MULTITFBS'} = $seventeen;
						$variantInformation{'ENSREGANNO_IMMUNOCELL_TFBS_SCORE'} = $eightteen;
						$variantInformation{'ENSREGANNO_IMMUNOCELL_MULTITFBS_SCORE'} = $nineteen;
					}
				
					my ($twenty, $twentyone) = PatternModule::RegAnno::getEvidenceType(\@per_cell_reg_features, $vf, "Multiple_Tissues");
					my ($twentytwo, $twentythree) = PatternModule::RegAnno::getEvidenceType(\@per_cell_reg_features, $vf, "Neuro");
					my ($twentyfour, $twentyfive) = PatternModule::RegAnno::getEvidenceType(\@per_cell_reg_features, $vf, "Immuno");
					$variantInformation{'ENSREGANNO_MULTICELL_EVIDENCE'} = $twenty;
					$variantInformation{'ENSREGANNO_MULTICELL_MULTIEVIDENCE'} = $twentyone;
					$variantInformation{'ENSREGANNO_NEUROCELL_EVIDENCE'} = $twentytwo;
					$variantInformation{'ENSREGANNO_NEUROCELL_MULTIEVIDENCE'} = $twentythree;
					$variantInformation{'ENSREGANNO_IMMUNOCELL_EVIDENCE'} = $twentyfour;
					$variantInformation{'ENSREGANNO_IMMUNOCELL_MULTIEVIDENCE'} = $twentyfive;
				}
			
				if ($variantInformation{'ENSREGANNO_REGULATORY_FEATURE'} eq "NA") {
					my ($twentysix, $twentyseven) = PatternModule::RegAnno::findClosestRegulatoryFeature($vf, $regfeat_adaptor, $slice_adaptor, $regint);
					$variantInformation{'ENSREGANNO_REGULATORY_FEATURE'} = $twentysix;
					$variantInformation{'ENSREGANNO_CONSEQUENCE'} = $twentyseven;
					if ($variantInformation{'ENSREGANNO_REGULATORY_FEATURE'} ne "NA") {
						my $reg_feature = $regfeat_adaptor->fetch_by_stable_id($variantInformation{'ENSREGANNO_REGULATORY_FEATURE'});
						my ($twentyeight, $twentynine) = PatternModule::RegAnno::findDistanceToStartAndEndOfRegFeature($reg_feature, $vf);
						$variantInformation{'ENSREGANNO_DISTANCE_TO_START'} = $twentyeight;
						$variantInformation{'ENSREGANNO_DISTANCE_TO_END'} = $twentynine;
					}
				}
			
				if ($variantInformation{'ENSREGANNO_REGULATORY_FEATURE'} ne "NA") {
					my @regFeat = split(" - ", $regfeat_adaptor->fetch_by_stable_id($variantInformation{'ENSREGANNO_REGULATORY_FEATURE'})->display_label());
					my @RegType = split(" ", $regFeat[0]);
					my $RegTypeNoWhite = join("_", @RegType);
					$variantInformation{'ENSREGANNO_REGULATORY_FEATURE_TYPE'} = $RegTypeNoWhite;
				}
			
				# If no overlapping TFBS site is found the closest TFBS is annotated, and hit is set to NO
				if ($variantInformation{'ENSREGANNO_MULTICELL_TFBS'} ne "NA") {
					$variantInformation{'ENSREGANNO_MULTICELL_TFBS_HIT'} = "YES";
				}
				if ($variantInformation{'ENSREGANNO_MULTICELL_TFBS'} eq "NA") {
					$variantInformation{'ENSREGANNO_MULTICELL_TFBS_HIT'} = "NO";
				}
				if ($variantInformation{'ENSREGANNO_MULTICELL_TFBS_HIT'} eq "NO") {
					$variantInformation{'ENSREGANNO_MULTICELL_TFBS'} = PatternModule::RegAnno::getNearestTranscriptionFactorBindingMotif($vf, $motifFeature_adaptor, $slice_adaptor, "Multiple_Tissues", $tfbsint);
				}
				if ($variantInformation{'ENSREGANNO_NEUROCELL_TFBS'} ne "NA") {
					$variantInformation{'ENSREGANNO_NEUROCELL_TFBS_HIT'} = "YES";
				}
				if ($variantInformation{'ENSREGANNO_NEUROCELL_TFBS'} eq "NA") {
					$variantInformation{'ENSREGANNO_NEUROCELL_TFBS_HIT'} = "NO";
				}
				if ($variantInformation{'ENSREGANNO_NEUROCELL_TFBS_HIT'} eq "NO") {
					$variantInformation{'ENSREGANNO_NEUROCELL_TFBS'} = PatternModule::RegAnno::getNearestTranscriptionFactorBindingMotif($vf, $motifFeature_adaptor, $slice_adaptor, "Neuro", $tfbsint);
				}
				if ($variantInformation{'ENSREGANNO_IMMUNOCELL_TFBS'} ne "NA") {
					$variantInformation{'ENSREGANNO_IMMUNOCELL_TFBS_HIT'} = "YES";
				}
				if ($variantInformation{'ENSREGANNO_IMMUNOCELL_TFBS'} eq "NA") {
					$variantInformation{'ENSREGANNO_IMMUNOCELL_TFBS_HIT'} = "NO";
				}
				if ($variantInformation{'ENSREGANNO_NEUROCELL_TFBS_HIT'} eq "NO") {
					$variantInformation{'ENSREGANNO_NEUROCELL_TFBS'} = PatternModule::RegAnno::getNearestTranscriptionFactorBindingMotif($vf, $motifFeature_adaptor, $slice_adaptor, "Immuno", $tfbsint);
				}
						
				# Printing annotation information to terminal
				print STDOUT $variantInformation{'ENSREGANNO_VARIANT_OVERLAPS_WITH_REGULATORY'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSREGANNO_REGULATORY_FEATURE'};
				print STDOUT "\t";		
				print STDOUT $variantInformation{'ENSREGANNO_REGULATORY_FEATURE_TYPE'};
				print STDOUT "\t";			
				print STDOUT $variantInformation{'ENSREGANNO_DISTANCE_TO_START'};
				print STDOUT "\t";		
				print STDOUT $variantInformation{'ENSREGANNO_DISTANCE_TO_END'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSREGANNO_CONSEQUENCE'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSREGANNO_MULTICELL_TFBS_HIT'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSREGANNO_MULTICELL_TFBS'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSREGANNO_MULTICELL_MULTITFBS'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSREGANNO_MULTICELL_TFBS_SCORE'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSREGANNO_MULTICELL_MULTITFBS_SCORE'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSREGANNO_MULTICELL_EVIDENCE'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSREGANNO_MULTICELL_MULTIEVIDENCE'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSREGANNO_NEUROCELL_TFBS_HIT'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSREGANNO_NEUROCELL_TFBS'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSREGANNO_NEUROCELL_MULTITFBS'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSREGANNO_NEUROCELL_TFBS_SCORE'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSREGANNO_NEUROCELL_MULTITFBS_SCORE'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSREGANNO_NEUROCELL_EVIDENCE'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSREGANNO_NEUROCELL_MULTIEVIDENCE'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSREGANNO_IMMUNOCELL_TFBS_HIT'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSREGANNO_IMMUNOCELL_TFBS'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSREGANNO_IMMUNOCELL_MULTITFBS'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSREGANNO_IMMUNOCELL_TFBS_SCORE'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSREGANNO_IMMUNOCELL_MULTITFBS_SCORE'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSREGANNO_IMMUNOCELL_EVIDENCE'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSREGANNO_IMMUNOCELL_MULTIEVIDENCE'};
				print STDOUT "\t";
			
				# Adding information to the output vcf-file
				$annotationString .= 
					";ENSREGANNO_VARIANT_OVERLAPS_WITH_REGULATORY=" . $variantInformation{'ENSREGANNO_VARIANT_OVERLAPS_WITH_REGULATORY'} .
					";ENSREGANNO_REGULATORY_FEATURE=" . $variantInformation{'ENSREGANNO_REGULATORY_FEATURE'} .
					";ENSREGANNO_REGULATORY_FEATURE_TYPE=" . $variantInformation{'ENSREGANNO_REGULATORY_FEATURE_TYPE'} .
					";ENSREGANNO_DISTANCE_TO_START=" . $variantInformation{'ENSREGANNO_DISTANCE_TO_START'} .
					";ENSREGANNO_DISTANCE_TO_END=" . $variantInformation{'ENSREGANNO_DISTANCE_TO_END'} .
					";ENSREGANNO_CONSEQUENCE=" . $variantInformation{'ENSREGANNO_CONSEQUENCE'} .
					";ENSREGANNO_MULTICELL_TFBS_HIT=" . $variantInformation{'ENSREGANNO_MULTICELL_TFBS_HIT'} .
					";ENSREGANNO_MULTICELL_TFBS=" . $variantInformation{'ENSREGANNO_MULTICELL_TFBS'} .
					";ENSREGANNO_MULTICELL_MULTITFBS=" . $variantInformation{'ENSREGANNO_MULTICELL_MULTITFBS'} .
					";ENSREGANNO_MULTICELL_TFBS_SCORE=" . $variantInformation{'ENSREGANNO_MULTICELL_TFBS_SCORE'} .
					";ENSREGANNO_MULTICELL_MULTITFBS_SCORE=" . $variantInformation{'ENSREGANNO_MULTICELL_MULTITFBS_SCORE'} .
					";ENSREGANNO_MULTICELL_EVIDENCE=" . $variantInformation{'ENSREGANNO_MULTICELL_EVIDENCE'} .
					";ENSREGANNO_MULTICELL_MULTIEVIDENCE=" . $variantInformation{'ENSREGANNO_MULTICELL_MULTIEVIDENCE'} .
					";ENSREGANNO_NEUROCELL_TFBS_HIT=" . $variantInformation{'ENSREGANNO_NEUROCELL_TFBS_HIT'} .
					";ENSREGANNO_NEUROCELL_TFBS=" . $variantInformation{'ENSREGANNO_NEUROCELL_TFBS'} .
					";ENSREGANNO_NEUROCELL_MULTITFBS=" . $variantInformation{'ENSREGANNO_NEUROCELL_MULTITFBS'} .
					";ENSREGANNO_NEUROCELL_TFBS_SCORE=" . $variantInformation{'ENSREGANNO_NEUROCELL_TFBS_SCORE'} .
					";ENSREGANNO_NEUROCELL_MULTITFBS_SCORE=" . $variantInformation{'ENSREGANNO_NEUROCELL_MULTITFBS_SCORE'} .
					";ENSREGANNO_NEUROCELL_EVIDENCE=" . $variantInformation{'ENSREGANNO_NEUROCELL_EVIDENCE'} .
					";ENSREGANNO_NEUROCELL_MULTIEVIDENCE=" . $variantInformation{'ENSREGANNO_NEUROCELL_MULTIEVIDENCE'} .
					";ENSREGANNO_IMMUNOCELL_TFBS_HIT=" . $variantInformation{'ENSREGANNO_IMMUNOCELL_TFBS_HIT'} .
					";ENSREGANNO_IMMUNOCELL_TFBS=" . $variantInformation{'ENSREGANNO_IMMUNOCELL_TFBS'} .
					";ENSREGANNO_IMMUNOCELL_MULTITFBS=" . $variantInformation{'ENSREGANNO_IMMUNOCELL_MULTITFBS'} .
					";ENSREGANNO_IMMUNOCELL_TFBS_SCORE=" . $variantInformation{'ENSREGANNO_IMMUNOCELL_TFBS_SCORE'} .
					";ENSREGANNO_IMMUNOCELL_MULTITFBS_SCORE=" . $variantInformation{'ENSREGANNO_IMMUNOCELL_MULTITFBS_SCORE'} .
					";ENSREGANNO_IMMUNOCELL_EVIDENCE=" . $variantInformation{'ENSREGANNO_IMMUNOCELL_EVIDENCE'} .
					";ENSREGANNO_IMMUNOCELL_MULTIEVIDENCE=" . $variantInformation{'ENSREGANNO_IMMUNOCELL_MULTIEVIDENCE'} ;
			}
			
			############################################################
			##### Annotation of transcription factor binding sites #####
			############################################################
			
			if ($allanno == 1 || $tfbsanno == 1) {
				my %variantInformation = PatternModule::TFBSAnno::createVariantTFBSHash();
				if ($vf) {
					
					# For variants that are not multiallelic
					if ($col[4] !~ /,/) {
						# All TFBSs overlapping the variant are annotated with position within the TFBS and the predicted effect on TF-binding
						my ($TFBSstring) = PatternModule::TFBSAnno::fetchVariantTFBSinfo($vf, $motifFeature_adaptor, $slice_adaptor, $csscore_adaptor, $mlss);  
						my ($TFBSconvertedString) = PatternModule::TFBSAnno::stringToHashConverter($TFBSstring, \%variantInformation);
						%variantInformation = %{$TFBSconvertedString};
					}
					
					# For variants that are multiallelic
					if ($col[4] =~ /,/) {
						# The variationfeature is divided into several variationfeatures
						my @variationFeatureList = PatternModule::GeneAnno::splitMultiAllelicVariationFeature($vf, $vf->slice(), $variant_feature_adaptor);

						# All TFBSs overlapping the variants are annotated with position within the TFBS and the predicted effect on TF-binding
						my ($TFBSstring) = PatternModule::TFBSAnno::fetchMultiallelicVariantTFBSinfo($vf, \@variationFeatureList, $motifFeature_adaptor, $slice_adaptor, $csscore_adaptor, $mlss);
						my ($TFBSconvertedString) = PatternModule::TFBSAnno::stringToHashConverter($TFBSstring, \%variantInformation);
						%variantInformation = %{$TFBSconvertedString};  
					}
				}
			
				# Printing annotation information to terminal
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
					print STDOUT $variantInformation{$TFBS1string};
					print STDOUT "\t";
					print STDOUT $variantInformation{$TFBS2string};
					print STDOUT "\t";
					print STDOUT $variantInformation{$TFBS3string};
					print STDOUT "\t";
					print STDOUT $variantInformation{$TFBS4string};
					print STDOUT "\t";
					print STDOUT $variantInformation{$TFBS5string};
					print STDOUT "\t";
					print STDOUT $variantInformation{$TFBS6string};
					print STDOUT "\t";
					print STDOUT $variantInformation{$TFBS7string};
					print STDOUT "\t";
					print STDOUT $variantInformation{$TFBS8string};
					print STDOUT "\t";
					print STDOUT $variantInformation{$TFBS9string};
					print STDOUT "\t";
					print STDOUT $variantInformation{$TFBS10string};
					print STDOUT "\t";
		
					# Counts what TFBS we are at
					$counter = $counter + 1;
				}							

				# Adding information to the output vcf-file
				$counter = 1;
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
					$annotationString .= 
						";" . $TFBS1string . "=" . $variantInformation{$TFBS1string} .
						";" . $TFBS2string . "=" . $variantInformation{$TFBS2string} .
						";" . $TFBS3string . "=" . $variantInformation{$TFBS3string} .
						";" . $TFBS4string . "=" . $variantInformation{$TFBS4string} .
						";" . $TFBS5string . "=" . $variantInformation{$TFBS5string} .
						";" . $TFBS6string . "=" . $variantInformation{$TFBS6string} .
						";" . $TFBS7string . "=" . $variantInformation{$TFBS7string} .
						";" . $TFBS8string . "=" . $variantInformation{$TFBS8string} .
						";" . $TFBS9string . "=" . $variantInformation{$TFBS9string} .
						";" . $TFBS10string . "=" . $variantInformation{$TFBS10string} ;
		
					# Counts what TFBS we are at
					$counter = $counter + 1;
				}
			}
		
			#######################################################
			##### Annotation of PolyPhen and SIFT predictions #####
			#######################################################
			if ($allanno == 1 || $predictanno == 1) {				
				my %variantInformation = ();
				
				if ($vf) {
					# Predicts most severe SIFT and PolyPhen consequences and returns hash with information
					(%variantInformation) = PatternModule::PredictAnno::predictConsequences($vf);
				}
				
				# Printing annotation information to terminal
				print STDOUT $variantInformation{'ENSPREDICTANNO_POLYPHEN_TRANSCRIPT_NAME'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSPREDICTANNO_POLYPHEN_PREDICTION'};
				print STDOUT "\t";		
				print STDOUT $variantInformation{'ENSPREDICTANNO_POLYPHEN_SCORE'};
				print STDOUT "\t";			
				print STDOUT $variantInformation{'ENSPREDICTANNO_SIFT_TRANSCRIPT_NAME'};
				print STDOUT "\t";		
				print STDOUT $variantInformation{'ENSPREDICTANNO_SIFT_PREDICTION'};
				print STDOUT "\t";
				print STDOUT $variantInformation{'ENSPREDICTANNO_SIFT_SCORE'};
				
				# Adding information to the output vcf-file
				$annotationString .= 
					";ENSPREDICTANNO_POLYPHEN_TRANSCRIPT_NAME=" . $variantInformation{'ENSPREDICTANNO_POLYPHEN_TRANSCRIPT_NAME'} .
					";ENSPREDICTANNO_POLYPHEN_PREDICTION=" . $variantInformation{'ENSPREDICTANNO_POLYPHEN_PREDICTION'} .
					";ENSPREDICTANNO_POLYPHEN_SCORE=" . $variantInformation{'ENSPREDICTANNO_POLYPHEN_SCORE'} .
					";ENSPREDICTANNO_SIFT_TRANSCRIPT_NAME=" . $variantInformation{'ENSPREDICTANNO_SIFT_TRANSCRIPT_NAME'} .
					";ENSPREDICTANNO_SIFT_PREDICTION=" . $variantInformation{'ENSPREDICTANNO_SIFT_PREDICTION'} .
					";ENSPREDICTANNO_SIFT_SCORE=" . $variantInformation{'ENSPREDICTANNO_SIFT_SCORE'} ;
			}
		}	

		############################
		##### Printing to file #####
		############################
		
		# For normal mode
		if ($isfix == 0) {
			$col[7] .= $annotationString;
		}
		
		# For fix mode
		if ($isfix == 1) {
		
			# If line needs to be re-annotated the info-section if fixed here, else the 
			# line is printed as it is 
			if ($annotate == 1) {
				my @info = split(";", $col[7]);
				
				# Array with indicies with information that need to be substituted
				my @subarray = ();

				# Opposite alphabetical order as GATK reorders annotations after individual split
				if ($lncanno == 1) {
					# Start and end of information to be substituted
					my $indexStart = firstidx { $_ =~ /ENSLNCANNO_GENE/ } @info;
					my $indexEnd = firstidx { $_ =~ /ENSLNCANNO_TRANSCRIPT_LENGTH/ } @info;
					
					# Array with the index range for substitutions
					my @subarray_lnc = ($indexStart .. $indexEnd);
					
					# Range added to overall substitution array
					push @subarray, @subarray_lnc;
				}
				if ($pcanno == 1) {
					# Start and end of information to be substituted
					my $indexStart = firstidx { $_ =~ /ENSPCANNO_GENE/ } @info;
					my $indexEnd = firstidx { $_ =~ /ENSPCANNO_VARIANT_OVERLAPS_WITH_GENE/ } @info;
					
					# Array with the index range for substitutions
					my @subarray_pc = ($indexStart .. $indexEnd);
					
					# Range added to overall substitution array
					push @subarray, @subarray_pc;
				}
				if ($reganno == 1) {
					# Start and end of information to be substituted
					my $indexStart = firstidx { $_ =~ /ENSREGANNO_CONSEQUENCE/ } @info;
					my $indexEnd = firstidx { $_ =~ /ENSREGANNO_VARIANT_OVERLAPS_WITH_REGULATORY/ } @info;
					
					# Array with the index range for substitutions
					my @subarray_reg = ($indexStart .. $indexEnd);
					
					# Range added to overall substitution array
					push @subarray, @subarray_reg;
				}
				if ($sncanno == 1) {
					# Start and end of information to be substituted
					my $indexStart = firstidx { $_ =~ /ENSSNCANNO_GENE/ } @info;
					my $indexEnd = firstidx { $_ =~ /ENSSNCANNO_TRANSCRIPT_LENGTH/ } @info;
					
					# Array with the index range for substitutions
					my @subarray_snc = ($indexStart .. $indexEnd);
					
					# Range added to overall substitution array
					push @subarray, @subarray_snc;
				}
				if ($tfbsanno == 1) {
					# Start and end of information to be substituted
					my $indexStart = firstidx { $_ =~ /ENSTFBSANNO_TFBS10_ALT_SCORE/ } @info;
					my $indexEnd = firstidx { $_ =~ /ENSTFBSANNO_TFBS9_SIZE_OF_INDEL/ } @info;
					
					# Array with the index range for substitutions
					my @subarray_tfbs = ($indexStart .. $indexEnd);
					
					# Range added to overall substitution array
					push @subarray, @subarray_tfbs;
				}
				if ($predictanno == 1) {
					# Start and end of information to be substituted
					my $indexStart = firstidx { $_ =~ /ENSPREDICTANNO_POLYPHEN_PREDICTION/ } @info;
					my $indexEnd = firstidx { $_ =~ /ENSPREDICTANNO_SIFT_TRANSCRIPT_NAME/ } @info;
					
					# Array with the index range for substitutions
					my @subarray_predict = ($indexStart .. $indexEnd);
					
					# Range added to overall substitution array
					push @subarray, @subarray_predict;
				}	
				if ($allanno == 1) {
					# Start and end of information to be substituted
					my $indexStart = firstidx { $_ =~ /ENSLNCANNO_GENE/ } @info;
					my $indexEnd = firstidx { $_ =~ /ENSTFBSANNO_TFBS9_SIZE_OF_INDEL/ } @info;
					
					# Array with the index range for substitutions
					my @subarray_all = ($indexStart .. $indexEnd);
					
					# Range added to overall substitution array
					push @subarray, @subarray_all;
				}
								
				# Testing whether there is information to be fixed
				die "\n\n\n\n\nInformation to be fixed is not present\n\n\nPlease try again!\n" unless (scalar @subarray > 0);
			
				# Removing information to be substituted and adding new information
				my %subarrayhash = map { $_ => 1 } @subarray;
				my @to_print = ();

				for (my $i = 0; $i < scalar @info; $i++) {
					if (!exists($subarrayhash{$i})) {
						push @to_print, $info[$i];
					}
				}

				$col[7] = join(";", @to_print) . 
						  $annotationString;
			}
		}
		
		print OUTPUT join("\t", @col), "\n";
		if ($annotate == 1) {
			print STDOUT "\n";
		}
	}	
} 
