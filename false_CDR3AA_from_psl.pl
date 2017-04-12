#!/usr/bin/perl
use strict;
use warnings;

## Identify CDR3 AA sequences that could be falsely identified by MiXCR
## Then a downstream remover can be applied to the MiXCR-exported *.clones.txt results

my $file1 = shift
  or die(
"Usage:$0 <*.alignments.txt> <blat_output.psl (against hg19/hg38)> <gene(TRA|TRB|TRG|TRD)>\n"
  );
my $file2 = shift
  or die(
"Usage:$0 <*.alignments.txt> <blat_output.psl (against hg19/hg38)> <gene(TRA|TRB|TRG|TRD)>\n"
  );
my $strand = shift
  or die(
"Usage:$0 <*.alignments.txt> <blat_output.psl (against hg19/hg38)> <gene(TRA|TRB|TRG|TRD)>\n"
  );

my %black_list;
my %white_list;
open IN, $file2 or die($!);
while (<IN>) {
    chomp;
    next unless $_;
    my @a = split;
    next unless $a[0] =~ /^\d+$/;
    my $outtarget = 1;
    if ( $strand eq "TRA" ) {
        $outtarget = 0
          if $a[13] =~ /chr14_KI27.*_alt/
              or (
                      $a[13] eq "chr14"
                  and $a[16] >= 21621838    # TRA locus start point on hg38
                  and $a[15] <= 23014042
              );                            # TRA locus end posint on hg19
    }
    elsif ( $strand eq "TRB" ) {
        $outtarget = 0

          #        print "$a[13]\n"
          if $a[13] =~
              /chr7_KI27.*_alt/ # chr7_KI270803v1_alt is the one carring paralogs
              or (
                      $a[13] eq "chr7"
                  and $a[16] >= 142000747    # TRB locus start position on hg19
                  and $a[15] <= 142813399
              );                             # TRB locus end position on hg38
    }
    elsif ( $strand eq "TRG" ) {
        $outtarget = 0
          if $a[13] =~ /chr7_KI27.*_alt/
              or (
                      $a[13] eq "chr7"
                  and $a[16] >= 38253380     # TRGJ locus start position on hg38
                  and $a[15] <= 38407494
              );                             # TRBV locus end position on hg19
    }
    elsif ( $strand eq "TRD" ) {
        $outtarget = 0
          if $a[13] =~ /chr14_KI27.*_alt/
              or (
                      $a[13] eq "chr14"
                  and $a[16] >= 22096032     # TRDV locus start position on hg38
                  and $a[15] <= 22928148
              );                             # TRDJ locus end position on hg19
    }
    else {
        die(
"Unrecognizable Gene name, which should be one of TRA, TRB, TRG or TRD\n"
        );
    }
    my $query_matched_proportion = ( $a[12] - $a[11] + 1 ) / $a[10];
    if ($outtarget) {
        next
          if $query_matched_proportion < 0.8
        ; # Do not consider query sequences with <80% matched to genome in blacklist
        next if exists $white_list{ $a[9] };
        $black_list{ $a[9] } = $a[0];
    }
    else {
        if ( exists $black_list{ $a[9] } ) {
            next if $black_list{ $a[9] } > $a[0];
            delete $black_list{ $a[9] };
        }
        $white_list{ $a[9] } = $a[0];
    }    # Filtering conditions

    #	print "$a[9]\t$a[13]\t$a[15]\t$a[16]\t$outtarget\n";
}
close IN;

my ( %fc3aa, %tc3aa );
my $n        = 0;
my $total_rc = 0;
my @cdr3s;
open IN, $file1 or die($!);
while (<IN>) {
    chomp;
    $n++;
    my @a = split /\t/;
    if ( $a[0] =~ /,/ ) {
        $fc3aa{ $a[5] }++
          if exists $black_list{"$n.1:$a[5]"}
              or exists $black_list{"$n.2:$a[5]"};
        $tc3aa{ $a[5] }++
          if exists $white_list{"$n.1:$a[5]"}
              or exists $white_list{"$n.2:$a[5]"};
    }
    else {
        $fc3aa{ $a[5] }++ if exists $black_list{"$n.0:$a[5]"};
        $tc3aa{ $a[5] }++ if exists $white_list{"$n.0:$a[5]"};
    }
}
close IN;

foreach my $k ( sort keys %fc3aa ) {
    $tc3aa{$k} = 0 unless exists $tc3aa{$k};
    print join "\t", $k, $tc3aa{$k}, $fc3aa{$k}, $tc3aa{$k} / $fc3aa{$k};
    print "\n";
}
