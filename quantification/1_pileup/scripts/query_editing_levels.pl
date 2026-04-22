#!/bin/perl
############################################################
#perl script that queries editing level of known sites in a BAM file

use warnings;
use strict;
use Getopt::Long;
require "parse_pileup_query.pl"; #NEED PARSE PILEUP LIBRARY

# Default values
my $minbasequal = 20;  # MINIMUM BASE QUALITY SCORE
my $minmapqual = 250;  # MINIMUM READ MAPPING QUALITY SCORE
my $offset = 33;       # BASE QUALITY SCORE OFFSET - 33 FOR SANGER SCALE, 64 FOR ILLUMINA SCALE
my $bamfile;
my $outputfile;
my $genomepath;        # PATH TO REFERENCE GENOME
my $inputfile;         # PATH TO EDITING SITE LIST (can be .bed or .bed.gz)
my $sampath = "samtools"; # PATH TO THE SAMTOOLS EXECUTABLE
my $help;

# Parse command line options
GetOptions(
    'bam|b=s'        => \$bamfile,
    'output|o=s'     => \$outputfile,
    'genome|g=s'     => \$genomepath,
    'sites|s=s'      => \$inputfile,
    'minbasequal|q=i' => \$minbasequal,
    'minmapqual|m=i'  => \$minmapqual,
    'offset|f=i'      => \$offset,
    'samtools|t=s'    => \$sampath,
    'help|h'         => \$help
) or die "Error in command line arguments\n";

# Help message
if ($help || !$bamfile || !$outputfile || !$genomepath || !$inputfile) {
    print "Usage: perl script.pl -b <bam> -o <output> -g <genome> -s <sites> [options]\n\n";
    print "Required arguments:\n";
    print "  -b, --bam          Input BAM alignment file (must be indexed)\n";
    print "  -o, --output       Output file name\n";
    print "  -g, --genome       Path to reference genome FASTA file\n";
    print "  -s, --sites        Path to editing site list (BED format, can be .bed.gz)\n\n";
    print "Optional arguments:\n";
    print "  -q, --minbasequal  Minimum base quality score (default: 20)\n";
    print "  -m, --minmapqual   Minimum read mapping quality score (default: 250)\n";
    print "                     Use 255 for unique mapping with STAR\n";
    print "                     Use ≥1 for reads mapped to <10 loci\n";
    print "  -f, --offset       Base quality score offset (default: 33)\n";
    print "                     Use 33 for Sanger/Illumina 1.8+\n";
    print "                     Use 64 for older Illumina scale\n";
    print "  -t, --samtools     Path to samtools executable (default: samtools)\n";
    print "  -h, --help         Show this help message\n\n";
    print "Example:\n";
    print "  perl script.pl -b sample.bam -o sample.rnaediting_op \\\n";
    print "    -g /path/to/genome.fasta -s /path/to/editing_sites.bed.gz\n\n";
    print "  perl script.pl -b sample.bam -o sample.out \\\n";
    print "    -g genome.fasta -s sites.bed -q 25 -m 255\n";
    exit(0);
}

# Verify BAM file exists
die "Error: BAM file does not exist: $bamfile\n" unless -e $bamfile;

# Verify BAM index exists
unless (-e "$bamfile.bai") {
    die "Error: BAM index not found. Please create index with: samtools index $bamfile\n";
}

# Verify required files exist
die "Error: Reference genome not found: $genomepath\n" unless -e $genomepath;
die "Error: Editing site list not found: $inputfile\n" unless -e $inputfile;

print STDERR "Creating temporary BED file...\n";
my $bedtemp = join '', $outputfile, '.bed';

# Check if input file is gzipped
my $awk_cmd;
if ($inputfile =~ /\.gz$/) {
    print STDERR "Detected gzipped input file, decompressing...\n";
    $awk_cmd = "zcat $inputfile | awk '\$1!=\"chromosome\"{print \$1\"\t\"\$2\"\t\"\$3}' > $bedtemp";
} else {
    $awk_cmd = "awk '\$1!=\"chromosome\"{print \$1\"\t\"\$2\"\t\"\$3}' $inputfile > $bedtemp";
}
system($awk_cmd);

print STDERR "Running samtools mpileup...\n";
my $piletemp = join '', $outputfile, '.pileup';
my $mpileup_cmd = "$sampath mpileup -A -B -d 1000000 -q $minmapqual -Q $minbasequal -f $genomepath -l $bedtemp $bamfile > $piletemp";
print STDERR "Command: $mpileup_cmd\n";
system($mpileup_cmd);

print STDERR "Parsing pileup output...\n";
my %sitehash;
my %genehash;  # NEW: Store gene IDs
open (my $PILEUP, "<", $piletemp) or die "Cannot open pileup file: $!\n";
my $pileup_lines = 0;
while(<$PILEUP>) {
	chomp;
	my ($chr, $position, $refnuc, $coverage, $pile, $qual) = split;
	my $location = join '_', $chr, $position;
	my ($refnuccount, $acount, $tcount, $ccount, $gcount) = &parse_pileup($_, $minbasequal, $offset);# parse each line of pileup
	my $counts = join ',', $refnuccount, $ccount, $gcount;
	$sitehash{$location} = $counts;
	$pileup_lines++;
}
close $PILEUP;
print STDERR "  Processed $pileup_lines pileup lines\n";

print STDERR "Cleaning up temporary files...\n";
system("rm $bedtemp");
system("rm $piletemp");

print STDERR "Writing output...\n";

# Open input file (handle gzipped or regular)
my $INPUT;
if ($inputfile =~ /\.gz$/) {
    open($INPUT, "zcat $inputfile |") or die "error opening gzipped inputfile: $!\n";
} else {
    open($INPUT, "<", $inputfile) or die "error opening inputfile: $!\n";
}

open (my $OUTPUT, ">", $outputfile) or die "error opening outputfile: $!\n";
print $OUTPUT "#chrom\tposition\tgene_id\tcoverage\teditedreads\teditlevel\n";  # CHANGED: Added gene_id column

my $sites_processed = 0;
my $sites_with_coverage = 0;

while (<$INPUT>) { #READ IN LIST OF KNOWN EDITED SITES AND QUERY EDITING STATUS
	chomp;
	my @fields = split;
	next if ($fields[0] eq 'chromosome');
	my ($chr, $position) = ($fields[0], $fields[2]); # 3rd column is 1-based coordinates
	my $location = join '_', $chr, $position;
	my ($strand) = ($fields[5]);
	my $gene_id = $fields[6];  # NEW: Get gene ID from column 7 (0-indexed as 6)

	$sites_processed++;

	if ($sitehash{$location}) { #PRINT OUT RESULT
		my ($refcount, $ccount, $gcount) = split(/\,/,$sitehash{$location});
		my ($newcov, $newmismatch) = (0,0);
		if ($strand eq '+') {
			$newmismatch = $gcount;
		} else {
			$newmismatch = $ccount;
		}
		$newcov = $refcount + $newmismatch;
		if ($newcov) {		
			my $varfreq = 0;
			$varfreq = sprintf("%.3f", $newmismatch/$newcov);
			print $OUTPUT "$fields[0]\t$fields[1]\t$gene_id\t$newcov\t$newmismatch\t$varfreq\n";  # CHANGED: Added $gene_id
			$sites_with_coverage++;
		} else {
			print $OUTPUT "$fields[0]\t$fields[1]\t$gene_id\t0\t0\tN/A\n";  # CHANGED: Added $gene_id
		}
	} else {
		print $OUTPUT "$fields[0]\t$fields[1]\t$gene_id\t0\t0\tN/A\n";  # CHANGED: Added $gene_id
	}
}
close $INPUT;	
close $OUTPUT;

print STDERR "\nDone!\n";
print STDERR "Summary:\n";
print STDERR "  Total sites queried: $sites_processed\n";
print STDERR "  Sites with coverage: $sites_with_coverage\n";
print STDERR "  Output saved to: $outputfile\n";