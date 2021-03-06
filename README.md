# MixcrFilter
[MiXCR](https://github.com/milaboratory/mixcr) is a powerful computation program for TCR analysis. It can identify human TCR clones from raw sequences in a highly sensitive way. However, when the input sequences were not exclusively from TCR regions, some paralog sequences, such as those from immunoglobulin regions, could be mis-identified as TCR genes. Here we provide a bash script, named MixcrFilter, to filter out TCR clones which were falsely recognized by MiXCR from paralog sequences according to alignment results by [BLAT](http://genome.ucsc.edu/goldenPath/help/blatSpec.html#blatUsage). 

## Pre-required programs:
1. [MiXCR](https://github.com/milaboratory/mixcr)
2. [BLAT](http://genome.ucsc.edu/goldenPath/help/blatSpec.html#blatUsage)

## Input
The sample name and paths to the input files need to be customized by editing Lines 5-10 of MixcrFilter.sh. The required input files are:

1. Raw RNA-seq reads in FASTQ format: `*.R1.fq.gz`, `*.R2.fq.gz`
2. Reference genome file, perferable hg38, which has better annotation of TCR sequences: `hg38.fa`

## Output
A sub-directary with the files with the sample name as the prefix, and following suffix:

| Suffix            | Description             |
| ----------------- | ----------------------- |
| alignments.vdjca  | Output from MiXCR align |
| clones.clns       | Output from MiXCR assemble |
| clones.txt        | Raw clones from MiXCR exportClones |
| cloneseq.fa       | FASTA file of all nucleotide sequences of CDR3 region |
| cloneseq.psl      | BLAT results of \*.cloneseq.fa against reference genome |
| blat1.clones.txt  | Clones excluding the ones which CDR3-containing reads mapped to non-TCR regions |
| blat2.clones.txt  | Clones excluding the ones which CDR3 nt seq mapped to non-TCR regions |
| alignments.txt    | Alignment data from MiXCR exportAlignments |
| alignedreads.fa   | Original reads containing identified CDR3 sequences | 
| alignedreads.psl  | BLAT results of \*.alignedreads.fa against reference genome |
| falseCDR3AA       | A list of AA sequences of CDR3 from falsely mapped reads (Format: AA-seq, \# In-target reads, \# Out-target reads, In/Out Ratio) |
| flt.clones.txt    | Filtered clones |





