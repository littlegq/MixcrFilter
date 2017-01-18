#!/usr/bin/perl
use strict;
use warnings;

# Remove recurrent CDR3 sequences, assuming CDR3 sequences should be unique among samples
# Also remove CDR3 sequences with less than 6 amino acid

my $min_cdr3aa = 6;   # remove CDR3 sequences with less than 6 amino acid

my $file1 = shift
  or die("Usage: $0 <*.blat.clones.txt> recurrent_CDR3.txt\n");
my $file2 = shift
  or die("Usage: $0 <*.blat.clones.txt> recurrent_CDR3.txt\n");
my $max_recur_criteria = 8;

my %black_list;
open IN, $file2 or die($!);
while (<IN>) {
    chomp;
    my @a = split;
    $black_list{ $a[0] } = 1 if $a[1] > $max_recur_criteria;
}
close IN;

my $total_rc = 0;
my @cdr3s;
open IN, $file1 or die($!);
while (<IN>) {
    chomp;
    my @a = split /\t/;
    next if @a < 8;   # remove incomplete records;
    if ( $a[0] =~ /^\d+$/ ) {
        next if exists $black_list{ $a[2] };
	   next if length($a[7]) < $min_cdr3aa;
        $total_rc += $a[0];
        push @cdr3s, $_;
    }
    else {
        print join "\t", @a;
        print "\n";
    }
}
close IN;

foreach my $cdr3 (@cdr3s) {
    my @a = split /\t/, $cdr3;
    $a[1] = $a[0] / $total_rc;
    print join "\t", @a;
    print "\n";
}

