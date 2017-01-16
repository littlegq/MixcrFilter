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
		$name/$name.$chain.clones.clns $name/$name.$chain.clones.txt
	
	echo "[`date`] Preparing FASTA files of CDR3 nt sequences as BLAT input for $pfx"
	tail -n +2 $name/$name.$chain.clones.txt | awk '{a++;print ">"a"\n"$3}' > $name/$name.$chain.cloneseq.fa
	
	echo "[`date`] Align CDR3 sequences of $pfx against human reference genome"
	blat $ref $name/$name.$chain.cloneseq.fa -ooc=$blat11ooc $name/$name.$chain.cloneseq.psl
	
	echo "[`date`] Filtering CDR3 clones based on blat results for $pfx "
	# remove clones with CDR3 sequences aligned to genomics regions other than expected TCR loci or,
	# with 90% or more of the sequences aligned to a single TCR gene
	./psl_modifier_for_mixcr.pl $name/$name.$chain.clones.txt $name/$name.$chain.cloneseq.psl $chain \
		> $name/$name.$chain.blat.clones.txt
	
	echo "[`date`] Preparing FASTA files of CDR3-containing reads as BLAT input for $pfx"
	$mixcr exportAlignments -sequence -vHit -jHit -vAlignment -jAlignment -aaFeature CDR3 \
		$name/$name.$chain.alignments.vdjca \
		$name/$name.$chain.alignments.txt
	awk 'NF==6' $name/$name.$chain.alignments.txt > $name/$name.$chain.alignments.txt.1
	mv $name/$name.$chain.alignments.txt.1 $name/$name.$chain.alignments.txt
	cat $name/$name.$chain.alignments.txt |\
		perl -ane '
			$n++;
			if( $F[0] =~ /,/ ){
				 my @seqs = split /,/, $F[0]; 
				 print ">$n.1:$F[5]\n$seqs[0]\n>$n.2:$F[5]\n$seqs[1]\n";
			} else {
			 	print ">$n.0:$F[5]\n$F[0]\n";
			}' > $name/$name.$chain.alignedreads.fa

	echo "[`date`] Align CDR3-containing reads of $pfx against human reference genome"
	blat $ref $name/$name.$chain.alignedreads.fa -ooc=$blat11ooc $name/$name.$chain.alignedreads.psl
	echo "[`date`] Identify and remove sequences that were falsely recognized as CDR3 for $pfx"
	./false_CDR3AA_remover.pl $name/$name.$gene.rec_cdr3.clones.txt $name/$name.$gene.falseCDR3AA \
		> $name/$name.$gene.falseCDR3AA
	./false_CDR3AA_remover.pl $name/$name.$gene.rec_cdr3.clones.txt $name/$name.$gene.falseCDR3AA \
		> $name/$name.$gene.flt.clones.txt
done

# move sample data to a sub-dicrectory
if [ ! -d "$pfx" ]
then
	mkdir $pfx
fi
mv $pfx.* $pfx


