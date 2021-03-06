#!/usr/bin/env perl

=head1 LICENSE

EphaGen Workflow Software
Copyright (c) 2017 Moscow Institute Of Physics and Technology

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Author: Maxim Ivanov

=head1 SYNOPSIS

ephagen.pl --bam FILE --ref FILE --vcf FILE --out FILE --out_vcf FILE

EphaGen - a package to estimate NGS Dataset quality in terms of 
its ability to detect mutations of predefined spectrum.
Given VCF file with mutations of interest it will calculate
probabilities of all mutation sites to be homozygous for reference
allele and based on these probabilities estimate Dataset sensitivity
to detect defined mutations. Input VCF file should contain ONLY mutations
of interest, while positions and reference alleles of defined mutations
should be in concordance with reference genome file.

=head1 ARGUMENTS

=over 9

=item --ref FILE

Path to reference genome fasta (required)

=item --bam FILE

Path to sample BAM file (required). BAM file should be aligned 
to the reference specified with "--ref" option, sorted and indexed.

=item --vcf FILE

Path to VCF file with target mutations. Note that mutation positions
in VCF file should be in concordance with input reference genome file.
Required unless --vcf_ref used.

=item --vcf_ref <BRCA|CFTR>

Instead of preparing VCF file with target mutations manually,
you can use one of the predefined files. 
Specify BRCA for pathogenic BRCA1/BRCA2 mutations 
(as of Breast Cancer Information core, https://research.nhgri.nih.gov/bic/, May 2018).
Specify CFTR for pathogenic CFTR mutation
(as of CFTR2 database, https://www.cftr2.org/, December 2017).
Note, that both files use hg19 human genome version as reference
so your input BAM file should be aligned to the same reference
For hg38 version of CFTR and BRCA variants files see ./reference directory.

=item --out FILE

Path to output file containing result dataset sensitivity analysis results.
This is tab-separated file where first raw contains information on mean coverage
and sensitivity for the whole dataset. Further lines contain results of 
downsample analysis based on random sampling of reads from the whole datasets.
For each read fraction random sampling is performed several times. Mean sensitivity 
is written for each read fraction. Mean coverage calculation 
is carried only across defined mutation sites.

=item --out_vcf FILE

Path to output VCF file containing sensitivity analysis per each mutation site

=item --sd [OPTIONAL]

Perform downsample analysis (default FALSE)

=item --p [OPTIONAL]

Number of cores to use

=item --da [OPTIONAL]

Degree of approximation of output sensitivity as number of decimals (default 3)

=back

=cut

use strict;
use warnings;
use Bio::DB::Sam;		#Install
use Bio::Cigar;			#Install
use Try::Tiny;			#Install
use Data::Dumper;
use List::Util qw(sum max min shuffle pairs);
use Switch;			#Install
use Storable 'dclone';
use Cwd qw(getcwd abs_path);
use Getopt::Long;
use File::Spec;
use Pod::Usage;
use File::Basename;
use Time::HiRes qw(gettimeofday);

$| = 1;

#-------------------------------------------------------------------------------------------
#---    CONSTANTS
#-------------------------------------------------------------------------------------------

my $version = '1.3';
my $qscore_min  = 15;
my $qscore_averaging_range      = 2;
my $var_qual    = 20;
my $downsample_av_number = 5;
my $mode;
my $downsample_config = "80/100;60/100;40/100;20/100;10/100;5/100;1/100";
my %vcf_definition = (
        'BRCA' => "BRCA.hg19.vcf",
        'CFTR' => "CFTR.hg19.vcf"
        );

#-------------------------------------------------------------------------------------------
#---    BEGIN
#-------------------------------------------------------------------------------------------

sub getAbsPath(\$) {
        my ($dirRef) = @_;
        my @tmp=glob($$dirRef);
        return 1 if(scalar(@tmp) != 1);
        my $ret = Cwd::realpath($tmp[0]);
        return 1 if !$ret && !($ret = File::Spec->rel2abs($tmp[0]));
        $$dirRef = $ret;
        return 0;
        }

my $baseDir;
my $refVCFDir;
my $tempDir;
my $statR_germline;
my $statR_somatic;
BEGIN {
	my $thisDir=(File::Spec->splitpath($0))[1];
	$thisDir = './' if ($thisDir eq "");
	$baseDir=File::Spec->catdir($thisDir,File::Spec->updir());
	$statR_germline=File::Spec->catdir($thisDir,'statG.r');
	$statR_somatic=File::Spec->catdir($thisDir,'statS.r');
	$refVCFDir=File::Spec->catdir($thisDir,File::Spec->updir(),'reference');
	$tempDir=File::Spec->catdir($thisDir,File::Spec->updir(),'temp');
	}

if(getAbsPath($baseDir)) {
	die "Can't resolve path for ephagen install directory: '$baseDir'\nExit status 1\n";
}

my $scriptName=(File::Spec->splitpath($0))[2];
my $argCount=scalar(@ARGV);
my $cmdline = join(' ',$0,@ARGV);

sub usage() { pod2usage(-verbose => 1,
		-exitval => 2); }

#-------------------------------------------------------------------------------------------
#---    USER CONFIGURATION
#-------------------------------------------------------------------------------------------

my ($inputBam, $refFile, $refVCF, $outFile, $outVCF, $vcfREF, $skip_downsample, $no_cores, $approx_degree);
my $help;

GetOptions( "bam=s" => \$inputBam,
	"ref=s" => \$refFile,
	"vcf=s" => \$refVCF,
	"vcf_ref=s" => \$vcfREF,
	"out=s" => \$outFile,
	"out_vcf=s" => \$outVCF,
	"somatic" => \$mode,
	"skip_downsample|sd" => \$skip_downsample,
	"p=s" => \$no_cores,
	"da=s" => \$approx_degree,
	"help|h" => \$help) || usage();

if (defined $mode) {$mode = 'somatic'} else {$mode = 'germline'}
$skip_downsample	//= 0;
$no_cores		//= 1;
$approx_degree		//= 2;

if (int($no_cores) ne $no_cores) {
	die "Can't resolve argument 'number of cores' - must be positive integer\nExit status 1\n";
	}
if ($no_cores < 1) {
	die "Can't resolve argument 'number of cores' - must be positive integer\nExit status 1\n";
	}

usage() if($help);
usage() unless($argCount);

#-------------------------------------------------------------------------------------------
#---    VALIDATE INPUT CONDITIONS
#-------------------------------------------------------------------------------------------

sub checkFile($;$) {
	my $file = shift;
	return if(-f $file);
	my $label = shift;
	die "Can't find" . (defined($label) ? " $label" : "") . " file: '$file'\nExit status 1\n";
}

sub checkFileArg($$) {
	my ($file,$label) = @_;
	die "Must specify $label file\nExit status 1\n" unless(defined($file));
	checkFile($file,$label);
	}

sub makeAbsoluteFilePaths(\$) {
	my ($filePathRef) = @_;
	my ($v,$fileDir,$fileName) = File::Spec->splitpath($$filePathRef);
	if (getAbsPath($fileDir)) {
#		die "Can't resolve directory path for '$fileDir' from input file argument: '$$filePathRef'\nExit status 1\n";
		}
	$$filePathRef = File::Spec->catfile($fileDir,$fileName);
	}

$statR_germline = abs_path($statR_germline);
makeAbsoluteFilePaths($statR_germline);
checkFileArg($statR_germline,"statG.r");

$statR_somatic = abs_path($statR_somatic);
makeAbsoluteFilePaths($statR_somatic);
checkFileArg($statR_somatic,"statS.r");


$inputBam = abs_path($inputBam);
makeAbsoluteFilePaths($inputBam);
checkFileArg($inputBam,"input BAM");

$refFile = abs_path($refFile);
makeAbsoluteFilePaths($refFile);
checkFileArg($refFile,"reference FASTA");

$tempDir = abs_path($tempDir);
makeAbsoluteFilePaths($tempDir);
if (opendir(TEMP, $tempDir)) {
	closedir TEMP;
	} else {
	if (mkdir($tempDir)) {
		} else {
		die "Can't find or create TEMP directory at $tempDir\nExit status 1\n";
		}
	}

my $indexR  = int(rand(89999)) + 10000;
my $Rinput  = File::Spec->catfile($tempDir, "Rinput$indexR");
my $Routput = File::Spec->catfile($tempDir, "Routput$indexR");
my $Rlog = File::Spec->catfile($tempDir, "Rlog$indexR");

if (defined($vcfREF)) {
	if (defined $vcf_definition{$vcfREF}) {
		$refVCF = File::Spec->catfile($refVCFDir, $vcf_definition{$vcfREF});
		makeAbsoluteFilePaths($refVCF);
		checkFileArg($refVCF,"reference $vcfREF VCF file");
		} else {
		die "Can't resolve value for option --vcf_ref. It should be one of the following: ",join(", ", keys %vcf_definition),"\nExit status 1\n";
		}
	} else {
	$refVCF = abs_path($refVCF);
	makeAbsoluteFilePaths($refVCF);
	checkFileArg($refVCF,"reference VCF file");
	}

$outFile = abs_path($outFile);
makeAbsoluteFilePaths($outFile);
if (open(TEST, ">", $outFile)) {close TEST} else {die "Can't use file $outFile as output\nExit status 1\n"}

$outVCF = abs_path($outVCF);
makeAbsoluteFilePaths($outVCF);
if (open(TEST, ">", $outVCF)) {close TEST} else {die "Can't use file $outVCF as output\nExit status 1\n"}

#-------------------------------------------------------------------------------------------
#---    CHECK BAM AND REFERENCE FASTA INDEX FILE
#-------------------------------------------------------------------------------------------

sub checkBamIndex($) {
	my ($file) = @_;
	my $ifile = $file . ".bai";
	if(! -f $ifile) {
		$ifile = $file;
		$ifile =~ s/\.bam$/\.bai/;
		if(! -f $ifile) {
			die "Can't find index for BAM file '$file'\nExit status 1\n";
			}
		}
	}

checkBamIndex($inputBam);

sub checkFaIndex($) {
	my ($file) = @_;
	my $ifile = $file . ".fai";
	if(! -f $ifile) {
		die "Can't find index for fasta file '$file'\nExit status 1\n";
		}
	# check that fai file isn't improperly formatted (a la the GATK bundle NCBI 37 fai files)
	open(my $FH,"< $ifile") || die "Can't open fai file '$ifile'\nExit status 1\n";
	my $lineno = 1;
	while(<$FH>) {
		chomp;
		my @F=split();
		if(scalar(@F) != 5) {
		die "Unexpected format for line number '$lineno' of fasta index file: '$ifile'\n\tRe-running fasta indexing may fix the issue. To do so, run: \"samtools faidx $file\"\nExit status 1\n";
		}
		$lineno++;
		}
	close($FH);
	}

checkFaIndex($refFile);

sub chechAssemblyConcordance {
	my $sam		= shift;
	my $vcfFile	= shift;
	my $refFile	= shift;
	
	my @chromosomes = $sam->features (-type=>'chromosome');
	my $faidx_file = $refFile . ".fai";
	
	open (my $VCF, "<", $vcfFile) || die "Can't open $vcfFile\nExit status 1\n";
	
	my @target_chr;
	while (<$VCF>) {
		chomp;
		next if m!^#!;
		my @fields = split/\t/;
		my $flag = 0;
		map {
			my $ref = $_;
#			print "$fields[0]\t$fields[1]\n";
			if ($ref->{'seqid'} eq $fields[0]) {
				die "Mutation '$fields[2]' start position is out of contig size at input VCF file $vcfFile\nExit status 1\n" if $fields[1] > $ref->{'end'};
				die "Mutation '$fields[2]' start position is out of contig size at input VCF file $vcfFile\nExit status 1\n" if $fields[1] < $ref->{'start'};
				$flag = 1;
				}
			} @chromosomes;
		die "Contig '$fields[0]' from input $vcfFile can't be found in input BAM file header\nExit status 1\n" if $flag eq 0;
		push (@target_chr, $fields[0]) unless grep(/^$fields[0]$/, @target_chr);
		}
	
	close($VCF);
	
	open (my $FAIDX, "<", $faidx_file) || die "Can't open $faidx_file\nExit status 1\n";
	
	while (<$FAIDX>) {
		chomp;
		my @fields = split/\t/;
		next unless grep(/^$fields[0]$/, @target_chr);
		my $flag = 0;
		map {
			my $ref = $_;
			if ($fields[0] eq $ref->{'seqid'}) {
				if ($fields[1] eq $ref->{'end'}) {
					$flag = 1;
					} else {die "Contig '$fields[0]' sizes from input BAM header and input reference file does not match\nExit status 1\n"}
				}
			} @chromosomes;
		die "Contig '$fields[0]' from input $refFile can't be found in input BAM file header\nExit status 1\n" if $flag eq 0;
		}
	
	close($FAIDX);
	
	}

#-------------------------------------------------------------------------------------------
#---    CORE INVOKE
#-------------------------------------------------------------------------------------------

head($inputBam, $refFile, $refVCF, $outFile, $outVCF, $version, $cmdline);

#-------------------------------------------------------------------------------------------
#---    CORE SUBROUTINES              
#-------------------------------------------------------------------------------------------


sub scigar {
	my $cigar	= shift;
	my $n		= shift;
	my $count_cigar	= 0;
	while ($n > 0) {
		if ( $cigar =~ /^(\d+)(S|M|I|D)/ ) {
			my $number = $1 - 1;
			my $letter = $2;
			if ( ( $letter eq "S" ) or ( $letter eq "I" ) ) {
				++$count_cigar;
				}
			if ( $number eq "0" ) {
				$cigar =~ s/^\d+\D//;
				} else {
					$cigar =~ s/^\d+/$number/;
					}
			} elsif ( $cigar =~ /^(\d+)(D)/m ) {
				$cigar =~ s/^\d+D//;
				next;
				}
		--$n;
		}
	return $count_cigar;
	}

sub slimRefAlt {
	my $REF = shift;
	my $ALT = shift;
	my $complete = shift;
	$complete = 0 unless defined $complete;
	my $delta = 0;
	my @ref = split//, $REF;
	my @alt = split//, $ALT;
	while ($ref[0] eq $alt[0]) {
		$delta += 1;
		@ref = join("", @ref[1..(scalar @ref - 1)]);
		@alt = join("", @alt[1..(scalar @alt - 1)]);
		}
	while ($ref[scalar @ref - 1] eq $alt[scalar @alt - 1]) {
		@ref = join("", @ref[0..(scalar @ref - 2)]);
		@alt = join("", @alt[0..(scalar @alt - 2)]);
		}
	$REF = join("", @ref); $ALT = join("", @alt);
	if ($complete eq 0) {
		$REF = '' if length($REF) < 1;
		$ALT = '' if length($ALT) < 1;
		}
	return ($REF, $ALT, $delta);
	}

sub conditioned_d_m {
	my $qual_mas	= shift;	# = [[ref], [alt]]
	my $m_ref	= shift;
	my $probability = 1;
	switch ($m_ref) {
		case 0 {	# 0 reference alleles - homozygous for alternative allele
			map {$probability = $probability * $_} @{@{$qual_mas}[0]};
			map {$probability = $probability * (1-$_)} @{@{$qual_mas}[1]};
			}
		case 1 {	# 1 reference allele  - heterozygote
			map {$probability = $probability / 2} @{@{$qual_mas}[0]};
			map {$probability = $probability / 2} @{@{$qual_mas}[1]};
			}
		case 2 {	# 2 reference alleles - homozygous for reference allele
			map {$probability = $probability * (1-$_)} @{@{$qual_mas}[0]};
			map {$probability = $probability * $_} @{@{$qual_mas}[1]};
			}
		}
	return $probability;
	}

sub conditioned_m_d {
	my $qual_mas	= shift;	# = [[ref], [alt]]
	my $m_ref	= shift;	
	my $prior_m	= shift;	# = [0, 1, 2]
	my $probability = conditioned_d_m($qual_mas, $m_ref) * @{$prior_m}[$m_ref];
	my $prior_d = 0;
	for (my $i = 0; $i <= 2; $i++) {
		$prior_d += conditioned_d_m($qual_mas, $i) * @{$prior_m}[$i];
		}
	if (	(@{$prior_m}[0] eq 1) and
		(@{$prior_m}[1] eq 0) and
		(@{$prior_m}[2] eq 0)
		) {
		return 0 if $m_ref eq 0;
		return $probability if $m_ref eq 1;
		return (1-$probability) if $m_ref eq 2;
		}
	if (($probability eq 0)and($prior_d eq 0)) {
		return @{$prior_m}[$m_ref];
		} else {
		die "Division by zero\nExit status 1\n" if $prior_d eq 0;
		$probability = $probability / $prior_d;
		return $probability;
		}
	}

sub score {
	my $phred	= shift;
	my $score = 10 ** (-$phred/10);
	return $score;
	}

sub unscore {
	my $score       = shift;
	return int(-10*log($score)/log(10))
	}

sub load_vcf {
	my $file	= shift;
	my %mut;

	open (READ, "<$file");

	while (<READ>) {
		chomp;
		next if m!^#!;
		my @mas = split/\t/;
		my $REF = $mas[3];
		die "Malformed reference allele at input VCF file $file for variant '$mas[2]'\nExit status 1\n" if length($REF) < 1;
		my $pos = $mas[1];
		my $ALT = $mas[4];
		die "Malformed reference allele at input VCF file $file for variant '$mas[2]'\nExit status 1\n" if length($ALT) < 1;
		my @alt = split/,/, $ALT;
		my @info = split/;/, $mas[7];
		for (my $alt_i = 0; $alt_i < scalar @alt; $alt_i++) {
			my $delta = 0;
			my $position = $mas[1];
			my $alt_slim;
			($mas[3], $alt_slim, $delta) = slimRefAlt($mas[3], $alt[$alt_i]);
			$position += $delta;
			my $count = 0;
			foreach my $arg (@info) {
				if ($arg =~ /COUNT=(\S+)/) {
					my @count_mas = split/,/, $1;
					die "Can't collect COUNT for each alternative allele for mutation '$mas[2]' in input VCF file $file\nExit status 1\n" if ((scalar @count_mas) ne (scalar @alt));
					$count = $count_mas[$alt_i];
					die "COUNT field should contain only INTEGER non-zero values\nExit status 1\n" if (int($count) ne $count);
					die "COUNT field should contain only INTEGER non-zero values\nExit status 1\n" if ($count eq 0);
					}
				}
			die "Can't collect COUNT for allele $ALT for variant $mas[2] in input VCF file $file\nExit status 1\n" if ($count eq 0);
			die "Can't collect COUNT for allele $ALT for variant $mas[2] in input VCF file $file\nExit status 1\n" if (int($count) ne $count);
			$mut{$mas[0]} = [@{$mut{$mas[0]}}, [$position, $mas[3], $alt_slim, [$mas[2], $pos, $REF, $alt[$alt_i], '', $mas[0]], $count, {}]] if defined $mut{$mas[0]};
			$mut{$mas[0]} = [[$position, $mas[3], $alt_slim, [$mas[2], $pos, $REF, $alt[$alt_i], '', $mas[0]], $count, {}]] unless defined $mut{$mas[0]};
			}
		}
	
	close READ;
	return \%mut;
	}

#mutation_hash {'seq_id'->[<element>]}
#element structure:
#	<0-based start position>,
#	<0-based reference allele>,
#	<0-based alternative allele>,
#	[
#		<mutation name>,
#		<1-based start position>,
#		<1-based reference allele>,
#		<1-based alternative allele>,
#		'',
#		<contig name>
#	],
#	<allele count>,
#	{
#		'GQ' -> [0, 1, 2],
#		'DP' -> [[base qualities for ref reads], [base qualities for alt reads]],
#		'DQ' -> [0, 1, 2],
#		'error' -> []
#		'name' -> ''
#		'prior' -> numeric
#		'USSR' -> numeric, single site sensitivity by R
#		'coverR' -> numeric, as calcuated by R script
#	}

sub get_allele_count {
	my $mutation_hash	= shift;
	
	my $count = 0;
	
	foreach my $seq_id (keys %{$mutation_hash}) {
		for (my $i = 0 ; $i < scalar @{($mutation_hash)->{$seq_id}}; $i++) {
			my $arg = (($mutation_hash)->{$seq_id})->[$i];
			$count  += ($arg)->[4];
			}
		}
	return $count;

	}

sub get_prior {
	my $mutation_hash	= shift;
	my $seq_id_ref		= shift;
	my $iArg		= shift;
	my $count = 0;
	my $af = 0;
	foreach my $seq_id (keys %{$mutation_hash}) {
		for (my $i = 0 ; $i < scalar @{($mutation_hash)->{$seq_id}}; $i++) {
			my $arg	= (($mutation_hash)->{$seq_id})->[$i];
			$count	+= ($arg)->[4];
			$af	+= ($arg)->[4] if (($i eq $iArg) and ($seq_id eq $seq_id_ref));
			}
		}
	my $ref_hom	= (1 - $af/$count)*(1 - $af/$count);
	my $hetero	= 2*($af/$count)*(1 - $af/$count);
	my $alt_hom	= ($af/$count)*($af/$count);
	return ($alt_hom, $hetero, $ref_hom);
	}

sub pipeline {
	my $mutation_hash	= shift;
	my $sam			= shift;
	my $mutation_hash_r	= dclone $mutation_hash;
	my $name		= 0;
	foreach my $seq_id (keys %{$mutation_hash_r}) {
		for (my $iArg = 0; $iArg < scalar @{($mutation_hash_r)->{$seq_id}}; $iArg++) {
			++$name;
			my $arg = (($mutation_hash_r)->{$seq_id})->[$iArg];
			my $start = ($arg)->[0];
			my $end = ($arg)->[0];
			$end += max(length(($arg)->[1]), length(($arg)->[2]));
			my $altR = ($arg)->[2];
			$altR =~ s/-//g;
			my $segment = $sam->segment($seq_id,$start,$end);
			my @all_alignments = $segment->features;
			my @dp = ([], []); # seq errors for each refs and alts
			foreach my $alignment (@all_alignments) {
				#aps,ape,rps,rpe,qps,qpe - 0-based coordinates
				#aps/ape - start/end position of variant site in alignment
				#rps/rpe - start/end positions of variation site in reference
				#qps/qpe - start/end positions of variation site in query(read)
				#warn "Could not parse tag values from alignment. Probably installed Bio::DB:Sam module is obsolete (required 1:43). Check your module version with \$cpan -D Bio::DB::Sam.";
				if (defined($alignment->get_tag_values("SUPPLEMENTARY"))) {
					next if $alignment->get_tag_values("SUPPLEMENTARY") eq '1';
					}
				next if $alignment->get_tag_values("UNMAPPED") eq '1';
				next if $alignment->get_tag_values("NOT_PRIMARY") eq '1';
				my $qps; my $rps; my $aps; my $qpe; my $rpe; my $ape;
				my $opqs; my $opqe;
				my ($ref, $match, $query) = $alignment->padded_alignment;
				my $cigar = Bio::Cigar->new($alignment->cigar_str);
				$rps = (($arg)->[0] - $alignment->start);
				$rpe = (($arg)->[0] - $alignment->start) + length(($arg)->[1]);
				my $dec = 0;
				next if $rps < 1;
				next if ($rps) > ($rpe);
				my $ts = $ref; $ts =~ s/-//g;next if ($rpe) > length($ts);
				try {($qps, $opqs) = $cigar->rpos_to_qpos($rps);};
				next unless defined $qps;
				$qpe = $qps + length(($arg)->[2]);
				$qpe -= 0;
				next unless defined $qpe;
				next if $qpe < 0;
				next if $qps < 0;
				$aps = $rps + scigar($alignment->cigar_str, $rps);
				$ape = $aps + max(length(($arg)->[2]),length(($arg)->[1]));
				my $oref = substr($ref, $aps, $ape-$aps);
				my $oalt = substr($query, $aps, $ape-$aps);
				my $omat = substr($match, $aps, $ape-$aps);
				$oref =~ s/-//g; $oalt =~ s/-//g;
				my @scores = @{$alignment->qscore};
				my $qscore_start = max(0, $qps - $qscore_averaging_range);
				my $qscore_end   = min($qpe - 1 + $qscore_averaging_range, (scalar @{$alignment->qscore}) - 1);
				my $qscore = 0;
				map {$qscore = $qscore + score($_)} @scores[($qscore_start)..($qscore_end)];
				$qscore = $qscore / ($qscore_end - $qscore_start + 1);
				next if $qscore >= score($qscore_min);
				if (($oref eq (($arg)->[1])) and ($oalt eq (($arg)->[2]))) {
					if (defined((($arg)->[5])->{'error'})) {
						if (grep(/mut/,@{(($arg)->[5])->{'error'}})) {
							} else {
                                                        push (@{$dp[1]},$qscore);
							}
						} else {
						push (@{$dp[1]},$qscore);
						}
					} elsif (((length($oref) eq length(($arg)->[1]))and
						(length($oalt) eq length(($arg)->[1])))or
						($omat =~ /^\|*$/)) {
						push (@{$dp[0]},$qscore);
						}
#				if ((($oref eq (($arg)->[1])) and ($oalt eq (($arg)->[1]))) or ($omat =~ /^\|*$/)) {
#					push (@{$dp[0]},$qscore);
#					} elsif (($oref eq (($arg)->[1])) and ($oalt eq (($arg)->[2]))) {
#					if (defined((($arg)->[5])->{'error'})) {
#						if (grep(/mut/,@{(($arg)->[5])->{'error'}})) {} else {
#							push (@{$dp[1]},$qscore);
#							}
#						} else {
#						push (@{$dp[1]},$qscore);
#						}
#					}
				}
			$dp[0] = [shuffle(@{$dp[0]})];
			$dp[1] = [shuffle(@{$dp[1]})];
			(((($mutation_hash_r)->{$seq_id})->[$iArg])->[5])->{'DP'} = \@dp;
			(((($mutation_hash_r)->{$seq_id})->[$iArg])->[5])->{'name'} = "m$name";
			}
		}
	return $mutation_hash_r;
	}

sub calling_wrapper { #calling subroutine wrapped to remove any systemic errors
	my $mutation_hash	= shift;
	my $mutation_hash_c	= $mutation_hash;
	foreach my $seq_id (keys %{$mutation_hash_c}) {
		for (my $iArg = 0; $iArg < scalar @{($mutation_hash_c)->{$seq_id}}; $iArg++) {
			calling($mutation_hash_c, $seq_id, $iArg);
			my $arg = (($mutation_hash_c)->{$seq_id})->[$iArg];
			if (defined((($arg)->[5])->{'error'})) { #if DQ or GQ indicates at systemic error -> purge all alternative alleles and repeat calling
				if (grep(/mut/,@{(($arg)->[5])->{'error'}})) {
					((((($mutation_hash_c)->{$seq_id})->[$iArg])->[5])->{'DP'})->[1] = [];
					calling($mutation_hash_c, $seq_id, $iArg);
					}
				}
			}
		}
	}

sub calling {
	my $mutation_hash	= shift;
	my $seq_id		= shift;
	my $mut_arg		= shift;
	
	my $mutation_hash_c = $mutation_hash;
	my $dp	= (((($mutation_hash_c)->{$seq_id})->[$mut_arg])->[5])->{'DP'};
	if (defined((((($mutation_hash_c)->{$seq_id})->[$mut_arg])->[5])->{'error'})) {
		delete ((((($mutation_hash_c)->{$seq_id})->[$mut_arg])->[5])->{'error'})
		};
	if (defined((((($mutation_hash_c)->{$seq_id})->[$mut_arg])->[5])->{'GQ'})) {
		delete ((((($mutation_hash_c)->{$seq_id})->[$mut_arg])->[5])->{'GQ'})
		};
	if (defined((((($mutation_hash_c)->{$seq_id})->[$mut_arg])->[5])->{'DQ'})) {
		delete ((((($mutation_hash_c)->{$seq_id})->[$mut_arg])->[5])->{'DQ'})
		};
	(((($mutation_hash_c)->{$seq_id})->[$mut_arg])->[5])->{'error'} = [];
	my $arg = (($mutation_hash_c)->{$seq_id})->[$mut_arg];
	my @prior = get_prior($mutation_hash_c, $seq_id, $mut_arg);
	my $hwe = 0;
	map {$hwe += conditioned_m_d($dp,$_,\@prior)} (0 ,1, 2);
	push (@{(((($mutation_hash_c)->{$seq_id})->[$mut_arg])->[5])->{'error'}}, 'hwe') unless $hwe eq 1;
	my @GQ = (conditioned_m_d($dp,0,\@prior), conditioned_m_d($dp,1,\@prior), conditioned_m_d($dp,2,\@prior));
	my @DQ = (conditioned_d_m($dp,0), conditioned_d_m($dp,1), conditioned_d_m($dp,2));
	(((($mutation_hash_c)->{$seq_id})->[$mut_arg])->[5])->{'GQ'} = \@GQ;
	(((($mutation_hash_c)->{$seq_id})->[$mut_arg])->[5])->{'DQ'} = \@DQ;
	my $var_qual_prob = (1 - (10**(-$var_qual/10)));
	if (((max @GQ) ne $GQ[2]) and ((max @GQ) > $var_qual_prob)) {
		push @{(((($mutation_hash_c)->{$seq_id})->[$mut_arg])->[5])->{'error'}}, 'mutG';
		} else {}
	if ((max @DQ) ne $DQ[2]) {
		push @{(((($mutation_hash_c)->{$seq_id})->[$mut_arg])->[5])->{'error'}}, 'mutD';
		}
	}

sub average_coverage {
	my $mutation_hash	= shift;
	
	my $sum = 0;
	my $count = 0;
	foreach my $seq_id (keys %{$mutation_hash}) {
		for (my $iArg = 0; $iArg < scalar @{($mutation_hash)->{$seq_id}}; $iArg++) {
			my $dp  = (((($mutation_hash)->{$seq_id})->[$iArg])->[5])->{'DP'};
			$sum += (scalar @{(($dp)->[0])});
			$sum += (scalar @{(($dp)->[1])});
			++$count;
			}
		}
	return $sum/$count;
	}

sub downsample {
	my $mutation_hash	= shift;
	my $lower		= shift;
	my $higher		= shift;
	my $mutation_hash_c	= dclone $mutation_hash;
	die "Undefined Error 'downsample1'\nExit status 1\n" unless $higher > $lower;
	my @ref_mas;
	for (my $i = 1; $i <= $lower; $i++) {
		push (@ref_mas, 1);
		}
	for (my $i = $lower + 1; $i <= $higher; $i++) {
		push (@ref_mas, 0);
		}
	@ref_mas = shuffle(@ref_mas);
	
	foreach my $seq_id (keys %{$mutation_hash_c}) {
		for (my $iArg = 0; $iArg < scalar @{($mutation_hash_c)->{$seq_id}}; $iArg++) {
			my $dp  = (((($mutation_hash_c)->{$seq_id})->[$iArg])->[5])->{'DP'}; # = [[ref], [alt]]
			my @ref;
			my @alt;
			@ref_mas = shuffle(@ref_mas);
			my $j = 0;
			for (my $i = 0; $i < scalar @ref_mas; $i++) {
				last unless defined ((($dp)->[0])->[$j]);
				push (@ref, ((($dp)->[0])->[$j])) if $ref_mas[$i] eq '1';
				++$j;
				$i = 0 if ($i eq ((scalar @ref_mas) - 1));
				}
			@ref_mas = shuffle(@ref_mas);
			for (my $i = 0; $i < scalar @ref_mas; $i++) {
				last unless defined ((($dp)->[0])->[$j]);
				push (@ref, ((($dp)->[0])->[$j])) if $ref_mas[$i] eq '1';
				++$j;
				$i = 0 if ($i eq ((scalar @ref_mas) - 1));
				}
			(((($mutation_hash_c)->{$seq_id})->[$iArg])->[5])->{'DP'} = [\@ref, \@alt];
			}
		}
	
	calling_wrapper($mutation_hash_c);
	return $mutation_hash_c;
	}

sub get_sens {
	my $mutation_hash	= shift;
	my $sens = 1;
	
	foreach my $seq_id (keys %{$mutation_hash}) {
		foreach my $arg (@{($mutation_hash)->{$seq_id}}) {
			if ( @{(($arg)->[5])->{'error'}} ) {
				} else {
				$sens = $sens * ((($arg)->[5])->{'GQ'})->[2];
				}
			}
		}
	
	return $sens;
	}

sub log10 {
	my $n		= shift;
	return log($n)/log(10);
	}

sub format_af {
	my $f		= shift;
	my $degree	= shift;
	$degree		//= 2;
	my $i = 1;
	my $c = 0;
	while (1) {
		if (int($i*$f) > (10**($degree - 1))) {
			my $s = log10($i);
			my $nom = sprintf "%.${s}f",$f;
			return $nom;
			}
		$i = $i * 10;
		++$c;
		return 0 if $c > 100;
#		die "Too small number '$f'\nExit status 1\n" if $c > 100;
		}
	}

sub format_sens {
	my $sens = shift;
	if ($sens =~ /^(0?\.?[9]+[^9]{1}).*/) {
#		$sens = $1;
		$sens = format_af($sens, $approx_degree);
		} else {
		$sens = format_af($sens, $approx_degree);
		}
	return $sens;
	}

sub write_vcf {
	my $mutation_hash	= shift;
	my $version		= shift;
	my $command		= shift;
	my $referencePath	= shift;
	my $fileHandler		= shift;
	
	print $fileHandler "##fileformat=VCFv4.2\n",
	"##ephaGenVersion=$version\n",
	"##ephaGenCommand=$command\n",
	"##reference=file:$referencePath\n",
	"##ALT=<ID=*,Description=\"Represents allele(s) other than observed.\">\n",
	"##INFO=<ID=ERROR,Number=0,Type=Flag,Description=\"Indicates that error occured during pileup or calling process.\">\n",
	"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read depth.\">\n",
	"##INFO=<ID=AF,Number=R,Type=Float,Description=\"Allele frequency.\">\n",
	"##INFO=<ID=SSS,Number=R,Type=Float,Description=\"Single Site Sensitivity\">\n",
	"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
	my $count = get_allele_count($mutation_hash);
	my @whole;
	map {push @whole, @{$mutation_hash->{$_}}} keys %{$mutation_hash};
#	foreach my $arg (sort {((($a)->[5])->{'GQ'})->[2] <=> ((($b)->[5])->{'GQ'})->[2]} @whole) {
	foreach my $arg (sort {(((($b)->[5])->{'prior'}) - ((($b)->[5])->{'USSR'})) <=> (((($a)->[5])->{'prior'}) - ((($a)->[5])->{'USSR'}))} @whole) {
		my $seq_id = ($arg->[3])->[5];
		my $line;
		my @info;
		my $warn_line = "Error at mutation site $seq_id\t".($arg->[3])->[1]."\t".($arg->[3])->[0]."\t".($arg->[3])->[2]."\t".($arg->[3])->[3];
		if (defined(((($arg)->[5])->{'coverR'}))) {
			push @info, "DP=".(((($arg)->[5])->{'coverR'}));
			} elsif (defined(((($arg)->[5])->{'DP'}))) {
				push @info, "DP=".((scalar @{((($arg)->[5])->{'DP'})->[0]}) + (scalar @{((($arg)->[5])->{'DP'})->[1]}));
				} else {
				warn "$warn_line\n";
				push @info, "ERROR";
				}
		push @info, "AF=".format_af(((($arg)->[4])/$count), 2);
		if ((defined(((($arg)->[5])->{'error'})))and(@{(($arg)->[5])->{'error'}})) {
			if ( grep(/hwe/,@{(($arg)->[5])->{'error'}}) )  {
				} elsif ( grep(/mut/,@{(($arg)->[5])->{'error'}}) ) {
				}
			warn "$warn_line\n";
			unless (grep(/ERROR/),@info) {push @info, "ERROR"};
			$line = join("\t", ($seq_id, ($arg->[3])->[1], ($arg->[3])->[0], ($arg->[3])->[2], ($arg->[3])->[3], '.', '.', join(";", @info)));
			} else {
#			push @info, "SSS=".format_sens(((($arg)->[5])->{'GQ'})->[2]);
			my $SSS = ((($arg)->[5])->{'USSR'})/((($arg)->[5])->{'prior'});
			$SSS = 1 if $SSS > 1;
			push @info, "SSS=".format_sens($SSS);
			$line = join("\t", ($seq_id, ($arg->[3])->[1], ($arg->[3])->[0], ($arg->[3])->[2], ($arg->[3])->[3], '.', '.', join(";", @info)));
			}
		print $fileHandler "$line\n";
		}
	}

sub downsample_wrapper {
	my $mutation_hash	= shift;
	my $config_string	= shift;
	my $downsample_count	= shift;
	my $fileHandler		= shift;
	
	my $string = $config_string;
	$string =~ s/\/|;/ /g;
	my @config = pairs(split/\s/, $string);
	
	my $duration = scalar @config;
	my $duration_l = length ($duration);
	my $j = 0;
	print STDERR "Downsample analysis:  ";
	printf STDERR "%${duration_l}d",0;
	print STDERR "/";
	printf STDERR "%${duration_l}d",$duration;
	map {
		my $low = $_->[0]; my $high = $_->[1];
		my @sens;
		my @avg_cov;
		
		my $mut_tmp = downsample($mutation_hash, $low, $high);
#		push (@sens, get_sens($mut_tmp));
		my $sens_tmp;
		($sens_tmp, $mut_tmp) = getSensR($mut_tmp);
		push (@sens, $sens_tmp);
		push (@avg_cov, average_coverage($mut_tmp));
		
		print $fileHandler "$low/$high\t";
		print $fileHandler format_af(&average(\@avg_cov), 2),"\t";
		print $fileHandler format_sens(&average(\@sens)),"\n";
		
		++$j;
		print STDERR "\r";
		print STDERR "Downsample analysis:  ";
		printf STDERR "%${duration_l}d",$j;
		print STDERR "/";
		printf STDERR "%${duration_l}d",$duration;
		} @config;
	print STDERR "\n";
	}

sub downsample_config_check {
	my $string		= shift;
	$string =~ s/\/|;/ /g;
	my @config = pairs(split/\s/, $string);
	map {die "downsample config: $_->[0] > $_->[1]\nExit status 1\n" if ($_->[0] > $_->[1])} @config;
	}

sub checkName {
	my $mutationHash	= shift;
	my %pikachu;
	
	my @whole_m;
	map {push @whole_m, @{$mutationHash->{$_}}} keys %{$mutationHash};
	foreach my $arg (@whole_m) {
		unless (defined((($arg->[5])->{"name"}))) {
			return 1;
			}
		my $name = (($arg->[5])->{"name"});
		if (defined($pikachu{$name})) {
			return 1;
			} else {
			$pikachu{$name} = 1;
			}
		}
	
	return 0;
	}

sub loadPrior {
	my $mutation_hash       = shift;
	my $count = 0;
	my $af = 0;
	
	my @whole_m;
	map {push @whole_m, @{$mutation_hash->{$_}}} keys %{$mutation_hash};
	
	foreach my $arg (@whole_m) {
		$count += $arg->[4];
		}
	
	die "Total allele count should be positive\nExit status 1\n" unless $count > 0;
	foreach my $seq_id (keys %{$mutation_hash}) {
		for (my $i = 0 ; $i < scalar @{($mutation_hash)->{$seq_id}}; $i++) {
			(((($mutation_hash)->{$seq_id})->[$i])->[5])->{"prior"} = ((($mutation_hash)->{$seq_id})->[$i])->[4]/$count;
                        }
                }
	}

sub getCountSum {
	my $mutation_hash	= shift;
	my $count = 0;
	
	my @whole_m;
	map {push @whole_m, @{$mutation_hash->{$_}}} keys %{$mutation_hash};
	
	foreach my $arg (@whole_m) {
		$count += $arg->[4];
		}
	
	return $count;
	}

sub loadR {
	my $mutationHash	= shift;
	my $RfileHandle		= shift;
	
	my @whole_m;
	my $count = getCountSum($mutationHash);
	my $countExp = 0;
	
	map {push @whole_m, @{$mutationHash->{$_}}} keys %{$mutationHash};
	my $priorSum += 0;
	foreach my $arg (@whole_m) {
#		next unless defined (($arg->[5])->{"DP"});
#		print Dumper $arg; exit;
		my $prior = ($arg->[4])/$count;;
		my $name = (($arg->[5])->{"name"});
		my @ref;
		my @alt;
		if (defined((($arg->[5])->{"DP"}))) {
			@ref  = @{(($arg->[5])->{"DP"})->[0]};
			@alt  = @{(($arg->[5])->{"DP"})->[1]};
			map {$_ = unscore($_)} @ref;
			map {$_ = unscore($_)} @alt;
			@ref = sort {$a <=> $b} @ref;
			@alt = sort {$a <=> $b} @alt;
			} else {
			@ref = ();
			@alt = ();
			}
		if ((defined($name))and(defined($prior))) {
			print $RfileHandle "$name\tprior\t$prior\n";
			print $RfileHandle "$name\tref\t",join("\t",@ref),"\n";
			print $RfileHandle "$name\talt\t",join("\t",@alt),"\n";
			$countExp += $arg->[4];
			} else {
			die "Undefined Error 'loadR1'\nExit status 1\n";
			}
		}
	die "Undefined Error 'loadR2'\nExit status 1\n" if ($count ne $countExp);
	}

sub readR {
	my $mutation_hash	= shift;
	my $RoutputHandle	= shift;
	
	my $sens;
	my %SSS;
	my %coverR;
	
	my $priorSum = 0;
	while (<$RoutputHandle>) {
		my $line = $_;
		chomp $line;
#		print STDERR "$line\n";
		if ($line =~ /^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+$/) {
#			print STDERR "$1\t$2\t$3\n";
			$SSS{$1} = $2;
			$priorSum += $3;
			$coverR{$1} = $4;
			} elsif ($line =~ /^(\S+)\s+$/) {
			$sens = $1;
			}
		}
#	die "Undefined Error 'readR1'- $priorSum\nExit status 1\n" unless $priorSum ne 1;
	
	foreach my $seq_id (keys %{$mutation_hash}) {
		for (my $i = 0 ; $i < scalar @{($mutation_hash)->{$seq_id}}; $i++) {
			my $name = (((($mutation_hash)->{$seq_id})->[$i])->[5])->{"name"};
			die "Undefined Error 'readR2' $Rinput $Routput\nExit status 1\n" unless (defined($SSS{$name}));
			(((($mutation_hash)->{$seq_id})->[$i])->[5])->{"USSR"} = $SSS{$name};
			(((($mutation_hash)->{$seq_id})->[$i])->[5])->{"coverR"} = $coverR{$name};
			}
		}
	return $sens;
	}

sub getSensR {
	my $mutationHash	= shift;
	my $mutationHashR	= dclone $mutationHash;
	
	open (my $RinputHandle, ">", $Rinput) or die "Can't open file $Rinput for writing\nExit status 1\n";
	loadR($mutationHashR, $RinputHandle);
	close $RinputHandle;
	
	print STDERR "Usage mode: $mode\n";
	if ($mode eq 'germline') {
		`R --slave -f $statR_germline --args $Rinput $no_cores > $Routput 2> $Rlog`;
		} elsif ($mode eq 'somatic') {
		`R --slave -f $statR_somatic --args $Rinput $no_cores > $Routput 2> $Rlog`;
		} else {
		die "Could not parse mode option (somatic or germline model)\n"
		}
	
	open (my $RoutputHandle, "<", $Routput) or die "Can't open file $Routput for reading\nExit status 1\n";
	my $sensR = readR($mutationHashR, $RoutputHandle);
	close $RoutputHandle;
	
	return ($sensR, $mutationHashR);
	}

#head($inputBam, $refFile, $refVCF, $outFile, $outVCF, $version, $cmdline);
sub head {
	my $inputBam	= shift;
	my $refFile	= shift;
	my $refVCF	= shift;
	my $outFile	= shift;
	my $outVCF	= shift;
	my $version	= shift;
	my $command	= shift;
	
	open (my $file_output, ">", $outFile) or die "Can't open file $outFile for writing\nExit status 1\n";
	open (my $file_vcf, ">", $outVCF) or die "Can't open file $outVCF for writing\nExit status 1\n";

	my $mut = load_vcf($refVCF);
	die "Empty VCF file\nExit status 1\n" if ((scalar (keys %{$mut})) eq 0);
	my @whole; map {push @whole, @{$mut->{$_}}} keys %{$mut};
	die "Provide 2 or more target variants in VCF file\nExit status 1\n" if ((scalar @whole) < 2);
	
	my $sam = Bio::DB::Sam->new(
		-bam   => "$inputBam",
		-fasta => "$refFile",
		-expand_flags  => 1,
		);
	
	chechAssemblyConcordance($sam, $refVCF, $refFile);
	downsample_config_check($downsample_config);
	
	print STDERR "Input parameters check: SUCCESS\n";
	print STDERR "Loading data from input BAM file...\n";
	
	my ($s1, $usec1) = gettimeofday();					#TIME
	my $mut1 = pipeline($mut, $sam);
	loadPrior($mut1);
	my ($s2, $usec2) = gettimeofday();					#TIME
	print STDOUT "TIME=",$s2-$s1+int(($usec2-$usec1)/10000)/100,"\n";	#TIME
	
	if (checkName($mut1)) {
		die "Undefined Error 'head1'\nExit status 1\n"
		};

	print STDERR "Calculating sensitivity\n";
	my $sensR;
	($sensR, $mut1) = getSensR($mut1);
	my ($s3, $usec3) = gettimeofday();					#TIME
	print STDOUT "TIME=",$s3-$s2+int(($usec3-$usec2)/10000)/100,"\n";	#TIME
	
#	print STDERR "Calculating sensitivity...\n";
#	calling_wrapper($mut1);
	
	print $file_output "#READ FRACTION\tMEAN COVERAGE\tSENSITIVITY\n";
	print $file_output "100/100\t";
	print $file_output format_af(average_coverage($mut1), 2),"\t";
#	print $file_output format_sens(get_sens($mut1)),"\t";
	print $file_output format_sens($sensR),"\n";
	
	downsample_wrapper($mut1, $downsample_config, $downsample_av_number, $file_output) if $skip_downsample;
	print STDERR "Writing VCF file...\n";
	write_vcf($mut1, $version, $command, $refFile, $file_vcf);	
	close $file_output;
	close $file_vcf;

	my $tmpFileRemove = deleteFile($Rinput) || deleteFile($Routput) || deleteFile($Rlog);
	warn "Can't remove temporary files at $tempDir\n" if $tmpFileRemove;
	}

sub deleteFile {
	my $file = shift;
	if (open(DELETEFILE, "<$file")) {
		close DELETEFILE;
		my $command;
		$command = `rm $file 2>&1`;
		chomp $command;
		if (length($command) > 1) {
			if (open(DELETEFILE, "<$file")) {
				close DELETEFILE;
				return 1;
				}
			}
		return 0;
		} else {
		return 0;
		}
	}

sub average{
	my($data) = @_;
	if (not @$data) {
		die("Empty array\nExit status 1\n");
		}
	my $total = 0;
	foreach (@$data) {
		$total += $_;
		}
	my $average = $total / @$data;
	return $average;
	}

sub stdev{
	my($data) = @_;
	if(@$data == 1){
		return 0;
		}
	my $average = &average($data);
	my $sqtotal = 0;
	foreach(@$data) {
		$sqtotal += ($average-$_) ** 2;
		}
	my $std = ($sqtotal / (@$data-1)) ** 0.5;
	return $std;
	}





















