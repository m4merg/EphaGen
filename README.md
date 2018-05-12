# EphaGen
**EphaGen - a package to estimate NGS Dataset qualitity in terms of its ability to detect mutations of predefined spectrum**

REQUIREMENTS

0. Perl installed and added to the PATH
1. SAMtools (http://samtools.sourceforge.net/, preferred version - 1.5) installed and added to the PATH
2. The following Perl libraries should be installed: Bio::DB::Sam (version 1.43), Bio::Cigar, Try::Tiny, Switch
3. The fasta sequence (one file, unzipped; e.g. "hg19.fa") of the targeted genome should be downloaded from http://hgdownload.soe.ucsc.edu/downloads.html and indexed with 'samtools faidx' command
4. You need to have your data aligned to target genome (.bam files). Bam files should be sorted and indexed with samtools commands 'samtools sort' and 'samtools index'
5. You need to define known mutations (VCF files) which will be used to calculate dataset sensitivity. Note, that VCF file should contain COUNT field in INFO section, representing allele count:
	INFO=<ID=COUNT,Number=R,Type=String,Description="Allele count">
   Note that input VCF should be in concordance with genome build you've used for data alignment.
   You can use one of the reference files stored in the ./reference directory: for pathogenic BRCA1/2 mutation analysis or pathogenic CFTR mutation analysis

INSTALLATION

0. Download EPHAGEN.zip
1. Unzip files
2. Check requirements (Perl libraries)
   You may use cpan to install the necessary Perl modules:
   
	   sudo cpan install Bio::DB::Sam
	   sudo cpan install Bio::Cigar
	   sudo cpan install Try::Tiny
	   sudo cpan install Switch
3. Run demo to check the installation
   Enter demo directory, invoke run.sh and compare output with example output from the corresponding directory

RUN EPHAGEN

To run ephagen execute src/ephagen.pl script. You may find example command line in demo/run.sh file.

---------------------------------------------------------------------------------------------------------------------------

****************LIMITATIONS OF USAGE****************

		Please make sure that:
	-	EphaGen handle only Single Nucleotide Variation, Multiple Nucleotide Variations (up to 50b.p. long), short insertions and deletions (up to 50b.p. long). Large genomic rearrangements, CNV, exon deletions and insertions are not supported.
	-	Reference VCF files, stored in the /reference directory (for BRCA and CFTR pathogenic mutation analysis) are made based on databased BREAST CANCER INFORMATION CORE and CFTR2. Allele counts taken from these databases refer to general population and may be inappropriate for some populations, especially for minor populations. Moreover, allele frequency spectrum based on these allele counts can possess bias towards overrepresented variations.

---------------------------------------------------------------------------------------------------------------------------

HOW TO READ OUTPUT FILES

EphaGen will produce two output files: general sensitivity analysis results in tsv format and VCF with probabilities for each mutation from input VCF file.
The first line from general sensitivity analysis results contains mean coverage across all mutation sites from input VCF and . Each further line contain results on downsample analysis: based on selected read fraction random sampling of reads from input BAM file is performed and sensitivity analysis invoked - mean sensitivity and standard deviation are written for each read fraction.
Output VCF file copies the input VCF files, but INFO section is rewritten.










