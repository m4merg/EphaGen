#Demo data to test functionality of EphaGen:
#You can find example output of this this script in ./example_output directory

#-------------------------------------------------------------------------------------------
#---    Change Paths If Needed
#-------------------------------------------------------------------------------------------

DATADIR=.
OUTDIR=.
TOOLDIR=..

cd $DATADIR
perl $TOOLDIR/src/ephagen.pl --bam $DATADIR/data/demo_data.bam --ref $DATADIR/ref/brac.fasta --vcf $DATADIR/BRAC.vcf --da 2 --out demo_out.tsv --out_vcf demo_out.vcf
perl check_result.pl
