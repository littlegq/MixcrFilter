#!/usr/bin/perl
use strict;
use warnings;

my $file1 = shift
  or die(
"Usage:$0 <*.clones.txt> <blat_output.psl (against hg19/hg38)> <chain(TRA|TRB|TRG|TRD)>\n"
  );
my $file2 = shift
  or die(
"Usage:$0 <*.clones.txt> <blat_output.psl (against hg19/hg38)> <chain(TRA|TRB|TRG|TRD)>\n"
  );
my $chain = shift
  or die(
"Usage:$0 <*.clones.txt> <blat_output.psl (against hg19/hg38)> <chain(TRA|TRB|TRG|TRD)>\n"
  );

my %black_list;
open IN, $file2 or die($!);
while (<IN>) {
    chomp;
    next unless $_;
    my @a = split;
    next unless $a[0] =~ /^\d+$/;
    my $outtarget = 1;
    if ( $chain eq "TRA" ) {
        $outtarget = 0
          if $a[13] eq "chr14"
              and $a[15] >= 21621838     # TRA locus start point on hg38
              and $a[16] <= 23014042;    # TRA locus end point on hg19
    }
    elsif ( $chain eq "TRB" ) {
        $outtarget = 0
          if $a[13] eq "chr7"
              and $a[15] >= 142000747     # TRB locus start position on hg19
              and $a[16] <= 142813399;    # TRB locus end position on hg38
    }
    elsif ( $chain eq "TRG" ) {
        $outtarget = 0
          if $a[13] eq "chr7"
              and $a[15] >= 38253380      # TRGJ locus start position on hg38
              and $a[16] <= 38407494;     # TRBV locus end position on hg19
    }
    elsif ( $chain eq "TRD" ) {
        $outtarget = 0
          if $a[13] eq "chr14"
              and $a[15] >= 22096032      # TRDV locus start position on hg38
              and $a[16] <= 22928148;     # TRDJ locus end position on hg19
    }
    else {
        die(
"Unrecognizable Gene name, which should be one of TRA, TRB, TRG or TRD\n"
        );
    }
    my $matched_size       = $a[12] - $a[11];
    my $matched_proportion = $matched_size / $a[10];
    my $match_score        = $a[0] / $a[10];
    $black_list{ $a[9] } = 1
      if ( $matched_proportion >= 0.9 and $match_score >= 0.9 )
      or $outtarget;                      # Filtering conditions
}
close IN;

my $n        = 0;
my $total_rc = 0;
my @cdr3s;
open IN, $file1 or die($!);
while (<IN>) {
    chomp;
    my @a = split /\t/;
    if ( $a[0] =~ /^\d+$/ ) {
        $n++;
        next if exists $black_list{$n};
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

