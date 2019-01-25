use strict;
use warnings;

if (open(READ, "<demo_out.tsv")) {close READ} else {die "Could not perform analysis\nExit\nError status: 1\n"}
if (open(READ, "<demo_out.vcf")) {close READ} else {die "Could not perform analysis\nExit\nError status: 1\n"}

open (READ, "<demo_out.vcf");

while (<READ>) {
	chomp;
	next if m!^#!;
	my @mas = split/\t/;
	my @info = split/;/, $mas[7];
	for (my $i = 0; $i <= 2; $i++) {
		if ($i eq 0) {
			if ($info[$i]=~/DP=(\d+)/) {
				} else {die "Could not perform analysis\nExit\nError status: 1\n"}
			}
		if ($i eq 1) {
			if ($info[$i]=~/AF=(\d+\.?\d*)/) {
				} else {die "Could not perform analysis\nExit\nError status: 1\n"}
			}
		if ($i eq 2) {
			if ($info[$i]=~/SSS=(\d+\.?\d*)/) {
				} else {die "Could not perform analysis\nExit\nError status: 1\n"}
			}
		}
	}

close READ;

open (READ, "<example_output/example_out.tsv");
my $head = <READ>;
my $line = <READ>;chomp $line;
my @mas = split/\t/, $line;
my $template_sens = $mas[2];
close READ;

open (READ, "<demo_out.tsv");
$head = <READ>;
$line = <READ>;
die "Could not perform analysis\nExit\nError status: 1\n" unless defined $line;
chomp $line;
die "Could not perform analysis\nExit\nError status: 1\n" if length($line) < 1;
@mas = split/\t/, $line;
my $sens = $mas[2];
die "WARNING: Analysis results are incorrect\nExit\nError status: 0\n" if (abs($sens - $template_sens)) > 0.05;
close READ;

die "Successfuly finished test\nExit\nError status: 0\n";




