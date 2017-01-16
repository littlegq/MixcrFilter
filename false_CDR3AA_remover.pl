#!/usr/bin/perl
use strict;
use warnings;

my $file1 = shift
  or die(
    "$0 <*.clones.txt> <*.falseCDR3AA>(output from false_CDR3AA_from_psl.pl)\n"
  );
my $file2 = shift
  or die(
    "$0 <*.clones.txt> <*.falseCDR3AA>(output from false_CDR3AA_from_psl.pl)\n"
  );

my %fc3aa;
open IN, $file2 or die($!);
while (<IN>) {
    chomp;
	my @a = split;
    $fc3aa{ $a[0] } = 1 if $a[3] < 1;   
		# Number of black_list counts must larger than the 
		# white_list counts from false_CDR3AA_from_psl.pl
}
close IN;

my $total_rc = 0;
my @cdr3s;
open IN, $file1 or die($!);
while (<IN>) {
    chomp;
	if(/^Clone_count/){
		print "$_\n";
		next;
	}
    my @a = split /\t/;
    next if exists $fc3aa{ $a[7] };
	$total_rc += $a[0];
	push @cdr3s, $_;
}
close IN;

foreach my $cdr3 (@cdr3s) {
	my @a = split /\t/, $cdr3;
	$a[1] = $a[0] / $total_rc;
	print join "\t", @a;
	print "\n";
}

