#!/usr/bin/perl -w

use strict;
use warnings;
use List::MoreUtils qw(firstidx);

=head1 NAME

Multi Allelic Line Splitter

=cut

=head1 DESCRIPTION

This script is designed to split lines containing multiple alternative alleles, coming
from multi sample calls, into lines containing only a single variant if possible. The 
genotypes are appropriately modified to fit the new structure.

=cut

=head1 SYNOPSIS

=head3 USAGE

perl MultiAllelicSplit.pl vcf_file_input.vcf vcf_file_output.vcf

=head3 INPUT

The script takes as input a VCF file that contains multiple alternative alleles, and 
contains the following in the info-section: "AC", "AF", "MLEAF" and "MLEAC".		
						
=head3 OUTPUT

The output is a vcf-file where the variants that contain multiple alternative alleles have
been split into multiple lines, with one alternative allele per line for 
those individuals that are homozygous, or those who are heterozygous for the reference. For
individuals that are heterozygous for two alternative alleles, a line with multiple 
alternative alleles is created. The genotypes are altered to fit the appropriate lines
(e.g. 0/3 in a line with multiple alternative alleles will be transformed into 0/1 in 
the line containing only the third alternative allele, while in the remaining lines it will 
be annotated as missing information "./.").									

=head1 FEEDBACK

Esben Eickhardt - Email E<lt>esbeneickhardt@biomed.au.dkE<gt>

=head1 AUTHORS

Esben Eickhardt - Email E<lt>esbeneickhardt@biomed.au.dkE<gt>

Francesco Lescai - Email E<lt>lescai@biomed.au.dkE<gt>

=cut

# Opening file and creating output file
open FILE, $ARGV[0] or die $!;
open OUTPUT, "> $ARGV[1]";

# Global variable and info-line to be added to the info section of the vcf-file
my $NumberOfIndividuals = 0;
my $info_line =  '##INFO=<ID=IS_MULTI_ALLELIC,Number=1,Type=String,Description="YES if both alleles are different from the reference">';


# Going through the file line by line
while (<FILE>) {
    chomp;
    my $line = $_;
    
    # Splitting line by tab
    my @col = split("\t", $line);
    my @col_unedited = split("\t", $line);
    
    # Printing header
    if ($line =~ /^##/) {
		if ($line =~ /^##INFO=<ID=MLEAC/) {
            print OUTPUT $info_line, "\n";
        }
        print OUTPUT $line, "\n";
        next;
    }
    
    # Counting number of individuals
    if ($line =~ /^#CHROM/) {
		$NumberOfIndividuals = (scalar @col) - 9;
        print OUTPUT $line, "\n";
        next;
    }      
    
    # Going through variant lines
    if ($line !~ /^#/) {
    	# If column number five does NOT contain a komma the line is printed as is
		if ($col[4] !~ /,/) {
			$col[7] .= ";IS_MULTI_ALLELIC=NO";
			print OUTPUT join("\t", @col), "\n";
		}
    
		# If column number five contains a komma the multi-allelic variant is split
		if ($col[4] =~ /,/) {
			# Creates a list of multi allelic variants
			my @multiAll = split(",", $col[4]);
			
			# Creates a list with all information in the info-section
			my @info = split(";", $col[7]);
			
			# Creates a list with all individuals Genotype information
			my @individuals = @col[9 .. (9+$NumberOfIndividuals-1)];
			
			# Empty list for placing genotypes of individuals where both alleles differ from the reference
			my @multiAllelicGenotypes = ();
			
			# Information on whether the line is complex will be added to this list
			my @is_complex = ();
			
			# Creates lists with information on MLEAF, MLEAC, AF and AC
			my @multiAC = SpecificInfoLister(\@info, "AC");
			my @multiAF = SpecificInfoLister(\@info, "AF");
			my @multiMLEAF = SpecificInfoLister(\@info, "MLEAF");
			my @multiMLEAC = SpecificInfoLister(\@info, "MLEAC");
	
			# Prints lines with single alternative alleles using a counter to indicate 
			# which allele is being processed
			for(my $i = 0; $i < scalar @multiAll; $i++) {
				# The alternative allele 
				$col[4] = $multiAll[$i];
				
				# Info section corrections
				$info[$multiAC[1]] = join("", "AC=", $multiAC[0][$i]);
				$info[$multiAF[1]] = join("", "AF=", $multiAF[0][$i]);
				$info[$multiMLEAF[1]] = join("", "MLEAF=", $multiMLEAF[0][$i]);
				$info[$multiMLEAC[1]] = join("", "MLEAC=", $multiMLEAC[0][$i]);
				$col[7] = join(";", @info);
				
				# Genotype adjustments				
				my @correctedGenotypes = ();
				@multiAllelicGenotypes = ();
				@is_complex = ();
				
				# For each individual the genotype is corrected matching the counter to
				# genotype
				for my $ind (@individuals) {
					my $corInd = $ind;
					# Ignoring missing information
					if (substr($ind, 1-1, 1) eq ".") {
						push @correctedGenotypes, $corInd;
						push @multiAllelicGenotypes, "./.";
						next;
					}
					
					# Ignoring 0/0 variants
					if (substr($ind, 1-1, 1) == 0 && substr($ind, 3-1, 1) == 0) {
						push @correctedGenotypes, $corInd;
						push @multiAllelicGenotypes, "./.";
						next;
					}
					
					# Catching complex variants with two non-reference
					if (substr($ind, 1-1, 1) > 0 && substr($ind, 3-1, 1) > 0 && substr($ind, 1-1, 1) != substr($ind, 3-1, 1)) {
						substr($corInd, 1-1, 1) = ".";
						substr($corInd, 3-1, 1) = ".";
						push @correctedGenotypes, $corInd;
						push @multiAllelicGenotypes, $ind;	
						push @is_complex, "YES";
						next;
					}
					
					# Catching homozygotes
					if (substr($ind, 1-1, 1) == substr($ind, 3-1, 1)) {
						if (substr($ind, 1-1, 1) == $i+1) {
							substr($corInd, 1-1, 1) = 1;
							substr($corInd, 3-1, 1) = 1;
						}
						if (substr($ind, 1-1, 1) != $i+1) {
							substr($corInd, 1-1, 1) = ".";
							substr($corInd, 3-1, 1) = ".";
						}					 
						push @correctedGenotypes, $corInd;
						push @multiAllelicGenotypes, "./.";
						next;					
					}
					
					# Catching heterozygotes
					if (substr($ind, 1-1, 1) == 0 || substr($ind, 3-1, 1) == 0) {					
						if (substr($ind, 1-1, 1) == $i+1) {
							substr($corInd, 1-1, 1) = 1;							
						}
						if (substr($ind, 3-1, 1) == $i+1) {
							substr($corInd, 3-1, 1) = 1;						
						}
						if (substr($ind, 1-1, 1) != $i+1) {
							substr($corInd, 1-1, 1) = 0;
						}
						if (substr($ind, 3-1, 1) != $i+1) {
							substr($corInd, 3-1, 1) = 0;
						}
						push @correctedGenotypes, $corInd;	
						push @multiAllelicGenotypes, "./.";					
					}					
				}
				
				# Adding individual information to columns
				my $Missing = grep { $_ eq "./." } @correctedGenotypes;
				my $refref = grep { $_ eq "0/0" } @correctedGenotypes;
				
				if ($Missing + $refref != scalar @correctedGenotypes) {
					@col[9 .. (9+$NumberOfIndividuals-1)] = @correctedGenotypes;
					$col[7] .= ";IS_MULTI_ALLELIC=NO";
				
					# Printing all information
					print OUTPUT join("\t", @col), "\n";
				}
			}
			if (scalar @is_complex > 0) {
				@col_unedited[9 .. (9+$NumberOfIndividuals-1)] = @multiAllelicGenotypes;
				$col_unedited[7] .= ";IS_MULTI_ALLELIC=YES";
				print OUTPUT join("\t", @col_unedited), "\n";
			}
		}
	}
}

# Pulls out a list of choice from the info section and the index of its position
sub SpecificInfoLister{
   my @list = @_;
   my @infoList = @{$list[0]};
   my @specInfo = grep(/^$list[1]=/, @infoList);
   my $index = firstidx { $_ =~ /$list[1]/ } @infoList;
   $specInfo[0] =~ s/$list[1]=//;
   my @multiSpecInfo = split(",", $specInfo[0]);
   return (\@multiSpecInfo, $index);
}

# Closing files
close(FILE);
close(OUTPUT);