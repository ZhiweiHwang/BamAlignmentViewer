# BamAlignmentViewer
Alignment Viewer of the special sites in bam base on python

```
###############################################################################
BamAlignmentViewer
A program for plot the reads covered target position for BAM sequence alignment
 files.
 Author: JoeyHwong (joeyhwong@hotmail.com)
 Licence: MIT
 usage: plotreads [-h] [--version] --bam [bamfile] --region [region] [options]
 positional arguments:
   -b, --bam     [FILE]  bam files containing mapped reads
   -r, --region  [STR]   the target region (0 base), example: chr1:123456-123457

optional arguments:
	-h, --help          show this help message and exit
	
	--version           show program's version number and exit
	--prefix    [STR]   string to add to the start of output graph.
	--rg        [STR]   The read group you want to draw
	-f          [NUM]   only show reads which mapping quality > NUM
	-d                  remove duplication reads, cigar level dups
	--ref       [FILE]  hg19 reference fasta file path
	-R                  show with reference, depend on -ref option specified

###############################################################################
```


