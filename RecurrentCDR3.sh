#!/bin/sh

# Stat the frequencies of CDR3 sequences
# Must run in the same directory as MixcrFilter.sh

cat */*.blat.clones.txt |\
	awk '$1!~/^Clone/' |\
	perl -e '
		my %cl;
		while (<>){
			chomp;
			my @a = split;
			$cl{$a[2]}++;
		}
		foreach my $cdr3 (sort keys %cl){
			print join "\t", $cdr3, $cl{$cdr3};
			print "\n";
		}' |\
	awk '$2>1' > recurrent_CDR3.txt


