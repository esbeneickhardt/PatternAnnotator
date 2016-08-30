# PatternAnnotator
Perl script for functionally annotating SNVs, short INDELs and multi-allelic variants using the Ensembl API

# Requirements
You must have an installation of the ENSEMBL API version 75 with the four databases: homo_sapiens_core_75_37, 
homo_sapiens_funcgen_75_37, homo_sapiens_variation_75_37 and ensembl_compara_75.

It is recommended that the databases are installed locally to avoid disconnects. Instructions for installing 
above mentioned can be found on the Ensembl website.

VCF-files must have the format outputted by GATK (version 2.7-2), so variants with alternative alleles in the 
<...> format seen in the 1000 genomes VCF-files cannot be handle.

# Running the script
By running the command "perldoc PatternAnnotator.pl" instructions of how to run the script and what it does
is available.

# Running on local databases
If you choose to run the script using local databases then the following needs to be changed such that you
query the local databases instead of the online databases:

    # For using online databases
    my $registry = 'Bio::EnsEMBL::Registry';
    
    $registry->load_registry_from_db(
        -host => 'ensembldb.ensembl.org',
        -user => 'anonymous'
    );

Instructions on how to query local databases through the API can also be found on the Ensembl website.
