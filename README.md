# EphaGen
**EphaGen - a package to estimate NGS Dataset qualitity in terms of its ability to detect mutations of predefined spectrum**

REQUIREMENTS

0. Perl installed and added to the PATH
1. SAMtools (http://samtools.sourceforge.net/, preferred version - 1.5) installed and added to the PATH
2. The following Perl libraries should be installed: Bio::DB::Sam (version 1.43), Bio::Cigar, Try::Tiny, Switch. The following R libraries should be installed: hash
3. The fasta sequence (one file, unzipped; e.g. "hg19.fa") of the targeted genome should be downloaded from http://hgdownload.soe.ucsc.edu/downloads.html and indexed with 'samtools faidx' command
4. You need to have your data aligned to target genome (.bam files). Bam files should be sorted and indexed with samtools commands 'samtools sort' and 'samtools index'
5. You need to define known mutations (VCF file) which will be used to calculate dataset sensitivity. Note, that VCF file should contain COUNT field in INFO section, representing allele count:
	INFO=<ID=COUNT,Number=R,Type=String,Description="Allele count">
   See example VCF files in the ./reference directory. Note that input VCF should be in concordance with genome build you've used for data alignment.
   You can use one of the reference files stored in the ./reference directory: for pathogenic BRCA1/2 mutation analysis or pathogenic CFTR mutation analysis

INSTALLATION

EphaGen is designed to be a stand alone command line application.

EphaGen requires the following system utilities: git, R, cpan. Check and make sure these are installed on your systems.

Go to a directory where EphaGen will be stored. Clone EphaGen from GitHub and then set EphaGen environment variable to point to where EphaGen is located

	   git clone https://github.com/m4merg/EphaGen
	   
	   cd EphaGen
	   
	   EPHAGEN_HOME=`pwd`

Install required perl and R libraries Check requirements:

	   wget http://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2
	   
	   tar xjf samtools-0.1.19.tar.bz2 && cd samtools-0.1.19 && make CFLAGS=-fPIC
	   
	   export SAMTOOLS=`pwd`
	   
	   sudo cpan install Bio::DB::Sam
	   
	   cd ../
	   
	   rm -r $SAMTOOLS
	   
	   sudo cpan install Bio::Cigar
	   
	   sudo cpan install Try::Tiny
	   
	   sudo cpan install Switch
	   
	   sudo su - -c "R -e \"install.packages('hash', repos='http://cran.rstudio.com/')\""

Run demo test to check the installation

	   cd $EPHAGEN_HOME/demo
	   
	   ./run.sh

Demo test should return Exit error status 0

RUN EPHAGEN

To run EphaGen execute:

	   perl $EPHAGEN_HOME/src/ephagen.pl

This will produce help message, describing required parameters and input files for running EphaGen.

Example usage:

	   ephagen.pl --bam FILE --ref FILE --vcf FILE --out FILE --out_vcf FILE
	   
Where:

	   --ref FILE
	   	Path to reference genome fasta (required)
		
	   --bam FILE
		Path to sample BAM file (required). BAM file should be aligned
		to the reference specified with "--ref" option, sorted and indexed.
		
	   --vcf FILE
		Path to VCF file with target mutations. Note that mutation positions
		in VCF file should be in concordance with input reference genome file.
		Required unless --vcf_ref used.
		
	   --out FILE
		Path to output file containing result dataset sensitivity analysis results.
		This is tab-separated file where first raw contains information on mean coverage
		and sensitivity for the whole dataset. Further lines contain results of
		downsample analysis based on random sampling of reads from the whole datasets.
		For each read fraction random sampling is performed several times. Mean sensitivity
		is written for each read fraction. Mean coverage calculation
		is carried only across defined mutation sites.
		
	   --out_vcf FILE
		Path to output VCF file containing sensitivity analysis per each mutation site

---------------------------------------------------------------------------------------------------------------------------

OUTPUT FILE DESCRIPTION

EphaGen will produce two output files: general sensitivity analysis results in tsv format defined by "--out" option and VCF with probabilities for each mutation from input VCF file defined by "--out_vcf" option.

General sensitivity analysis results provide sensitivity calculation results for each downsample iteration line by line. If no downsample was carried out, file contains only one line for 100% of reads. In addition to sensitivity mean coverage is provided.
Output VCF file copies the input VCF files provided with "--vcf" option with INFO section is expanded. Field "SSS" is added into INFO section providing information for Single Site Sensitivity for each mutation:

	   ##INFO=<ID=SSS,Number=R,Type=Float,Description="Single Site Sensitivity">
	   
Variants in the output VCF file are sorted based on the false negative rate in descending order.

---------------------------------------------------------------------------------------------------------------------------

****************LIMITATIONS OF USAGE****************

Please make sure that:
-	EphaGen handle only Single Nucleotide Variation, Multiple Nucleotide Variations (up to 50b.p. long), short insertions and deletions (up to 50b.p. long). Large genomic rearrangements, CNV, exon deletions and insertions are not supported.
-	Reference VCF files, stored in the /reference directory (for BRCA and CFTR pathogenic mutation analysis) are made based on databases BREAST CANCER INFORMATION CORE and CFTR2. Allele counts taken from these databases refer to general population and may be inappropriate for some populations, especially for minor populations. Moreover, allele frequency spectrum based on these allele counts can possess bias towards overrepresented variations.
