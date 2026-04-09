use strict;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use File::Path qw(make_path);
use File::Basename qw(dirname basename);
use Getopt::Long;

my %sitehash;
my %totalhash;
my %lvlhash;

# Default values
my $minsamps = 450;
my $mincov = 20;
my $input_dir;
my $output_file;
my $help;

# Parse command line options
GetOptions(
    'input|i=s'      => \$input_dir,
    'output|o=s'     => \$output_file,
    'mincov|c=i'     => \$mincov,
    'minsamps|s=i'   => \$minsamps,
    'help|h'         => \$help
) or die "Error in command line arguments\n";

# Help message
if ($help || !$input_dir || !$output_file) {
    print "Usage: perl script.pl -i <input_dir> -o <output_file> [options]\n\n";
    print "Required arguments:\n";
    print "  -i, --input     Input directory containing .rnaediting_op files\n";
    print "  -o, --output    Output file path\n\n";
    print "Optional arguments:\n";
    print "  -c, --mincov    Minimum coverage per site (default: 20)\n";
    print "  -s, --minsamps  Minimum number of samples per site (default: 450)\n";
    print "  -h, --help      Show this help message\n\n";
    print "Example:\n";
    print "  perl script.pl -i /path/to/input -o /path/to/output.txt\n";
    print "  perl script.pl -i /path/to/input -o /path/to/output.txt -c 10 -s 100\n";
    exit(0);
}

# Verify input directory exists
die "Error: Input directory does not exist: $input_dir\n" unless -d $input_dir;

# Create output directory if it doesn't exist
my $output_dir = dirname($output_file);
make_path($output_dir) unless -d $output_dir;

print STDERR "Parameters:\n";
print STDERR "  Input directory: $input_dir\n";
print STDERR "  Output file: $output_file\n";
print STDERR "  Minimum coverage: $mincov\n";
print STDERR "  Minimum samples: $minsamps\n\n";

open(my $OUT, ">", $output_file) or die "Cannot open output file: $!\n";

# Process each input file
my $file_count = 0;
foreach my $file (glob "$input_dir/*.rnaediting_op") {
    $file_count++;
    print STDERR "Processing file $file_count: $file\n";

    # sample = filename only (no path, no extension)
    my $sample = basename($file, ".rnaediting_op");
    print STDERR "  Sample name: $sample\n";

    open (my $INPUT, "<", $file);
    while(<$INPUT>) {
        chomp;
        my @fields = split;
        my ($chr, $pos, $cov, $edit, $lvl) = ($fields[0],$fields[1],$fields[2], $fields[3],$fields[4]);
        my $site = join ':', $chr, ($pos-1), $pos;
        my $ratio = join '/', $edit, $cov;
        next if ($chr eq '#chrom');
        if ($cov >= $mincov) {
            $sitehash{$sample}{$site} = $ratio;
            if ($totalhash{$site}) {
                $totalhash{$site}++;
                $lvlhash{$site} = join ',', $lvlhash{$site}, $ratio;
            } else {
                $totalhash{$site} = 1;
                $lvlhash{$site} = $ratio;
            }
        }
    }
    close $INPUT;
}

print STDERR "\nTotal samples processed: $file_count\n";
print STDERR "Writing output matrix...\n";

# Write header
print $OUT "chrom";
foreach my $sample (keys %sitehash) {
    print $OUT " $sample";
}
print $OUT "\n";

# Write data rows
my $sites_written = 0;
my $sites_filtered = 0;
foreach my $site (keys %totalhash) {
    if ($totalhash{$site} >= $minsamps) {
        my @lvls = split(/\,/, $lvlhash{$site});
        @lvls = sort {$b <=> $a} @lvls;
        print $OUT "$site";
        foreach my $sample (keys %sitehash) {
            if ($sitehash{$sample}{$site}) {
                print $OUT " $sitehash{$sample}{$site}";
            } else {
                print $OUT " 0/0";
            }
        }
        print $OUT "\n";
        $sites_written++;
    } else {
        $sites_filtered++;
    }
}

close $OUT;

print STDERR "\nDone!\n";
print STDERR "Sites written: $sites_written\n";
print STDERR "Sites filtered (< $minsamps samples): $sites_filtered\n";
print STDERR "Output saved to: $output_file\n";