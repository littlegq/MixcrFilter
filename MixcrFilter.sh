#!/bin/bash
# Author: Qiang Gong <qgong@coh.org>

# Before running, please customize/examine the following options:
pfx="SampleID"
RawDataDir="./Rawdata"
r1="$RawDataDir/$pfx.R1.fastq.gz"
r2="$RawDataDir/$pfx.R2.fastq.gz"
ref="hg38.fa"
blat11ooc="hg38.11.ooc"
if [ ! -f $blat11ooc ]
then
	blat $ref /dev/null /dev/null -makeOoc=$blat11ooc -repMatch=1024
fi

######################### 


for chain in "TRA" "TRB" "TRG" "TRD"
do
	echo "[`date`] Call $chain of $pfx"
	
	echo "[`date`] Building alignments using mixcr for $pfx"
	mixcr align --loci $chain $r1 $r2 $pfx.$chain.alignments.vdjca
	
	echo "[`date`] Assembling clones using mixcr for $pfx "
	mixcr assemble $pfx.$chain.alignments.vdjca $pfx.$chain.clones.clns
	
	echo "[`date`] Exporting clones to tab-delimited files using mixcr for $pfx"
	mixcr exportClones -count -fraction -sequence -vHit -jHit -vAlignment -jAlignment -aaFeature CDR3 \
		$pfx.$chain.clones.clns $pfx.$chain.clones.txt
	
	echo "[`date`] Preparing FASTA files of CDR3 nt sequences as BLAT input for $pfx"
	tail -n +2 $pfx.$chain.clones.txt | awk '{a++;print ">"a"\n"$3}' > $pfx.$chain.cloneseq.fa
	
	echo "[`date`] Align CDR3 sequences of $pfx against human reference genome"
	blat $ref $pfx.$chain.cloneseq.fa -ooc=$blat11ooc $pfx.$chain.cloneseq.psl
	
	echo "[`date`] Filtering CDR3 clones based on blat results for $pfx "
	# remove clones with CDR3 sequences aligned to genomics regions other than expected TCR loci or,
	# with 90% or more of the sequences aligned to a single TCR gene
	./psl_modifier_for_mixcr.pl $pfx.$chain.clones.txt $pfx.$chain.cloneseq.psl $chain \
		> $pfx.$chain.blat.clones.txt
	
	echo "[`date`] Preparing FASTA files of CDR3-containing reads as BLAT input for $pfx"
	$mixcr exportAlignments -sequence -vHit -jHit -vAlignment -jAlignment -aaFeature CDR3 \
		$pfx.$chain.alignments.vdjca \
		$pfx.$chain.alignments.txt
	awk 'NF==6' $pfx.$chain.alignments.txt > $pfx.$chain.alignments.txt.1
	mv $pfx.$chain.alignments.txt.1 $pfx.$chain.alignments.txt
	cat $pfx.$chain.alignments.txt |\
		perl -ane '
			$n++;
			if( $F[0] =~ /,/ ){
				 my @seqs = split /,/, $F[0]; 
				 print ">$n.1:$F[5]\n$seqs[0]\n>$n.2:$F[5]\n$seqs[1]\n";
			} else {
			 	print ">$n.0:$F[5]\n$F[0]\n";
			}' > $pfx.$chain.alignedreads.fa
	
	echo "[`date`] Align CDR3-containing reads of $pfx against human reference genome"
	blat $ref $pfx.$chain.alignedreads.fa -ooc=$blat11ooc $pfx.$chain.alignedreads.psl
	./false_CDR3AA_from_psl.pl $pfx.$chain.blat.clones.txt $pfx.$chain.falseCDR3AA \
		> $pfx.$chain.falseCDR3AA
	./false_CDR3AA_remover.pl $pfx.$chain.blat.clones.txt $pfx.$chain.falseCDR3AA \
		> $pfx.$chain.blat1.clones.txt
	tail -n +2 $pfx.$chain.blat1.clones.txt | awk '{a++;print ">"a"\n"$3}' > $pfx.$chain.cloneseq.fa
	blat $ref $pfx.$chain.cloneseq.fa -ooc=$blat11ooc $pfx.$chain.cloneseq.psl
	./psl_modifier_for_mixcr.pl $pfx.$chain.blat1.clones.txt $pfx.$chain.cloneseq.psl $chain > $pfx.$chain.blat2.clones.txt
	
	## *blat2.clones.txt files from different samples could be used to summarize the frequencies of CDR3 sequences (see RecurrentCDR3.sh)
	## We provide our observed frequency in the file recurrent_CDR3.txt
	if [ -f recurrent_CDR3.txt ]
	then
		echo "[`date`] Removing clones with recurrent CDR3 sequences for $pfx"
		./recurrent_removal.pl $pfx.$chain.blat2.clones.txt recurrent_CDR3.txt |\
			sed 's/ /_/g' > $pfx.$chain.flt.clones.txt
	else
		cp $pfx.$chain.blat2.clones.txt $pfx.$chain.flt.clones.txt
	fi
done

# move sample data to a sub-dicrectory
if [ ! -d "$pfx" ]
then
	mkdir $pfx
fi
mv $pfx.* $pfx


