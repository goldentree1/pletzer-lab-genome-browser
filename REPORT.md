# Elliott Brown Report

My time at Pletzer Lab has been an interesting journey! 
Me, a computer scientist, working in the Microbiology building for a Summer Student Internship job.
The job had been advertised as creating interactive website charts for data, which I have had an interest in, having made a surf forecast project of my own (of a similar nature, though this was certainly a much bigger challenge).

At first - lacking a lot of fundamental knowledge of genetics (and biology in general) - my first assignment was to learn. I became accustomed to the base nucleotides (ACGT/U) that make up DNA and RNA, how RNA is copied from DNA and can be used as messenger for gene expression, how gene expression can give different experimental results, for example, an infected culture may express different genes to a control group. It was fascinating learning about the consequences of this: how we can use these results to study bacterial resistance to antibiotics.

More relevant to my work was the raw data associated with RNA-sequencing (RNA-seq), and how this was processed. 
Many such files exist:
  - FASTA/FASTQ files, which are simple and just describe a sequence of nucleotides at certain read positions. These were used in my project as the reference genome - the DNA sequence that is well known, studied and confirmed in the NCBI database. We can use a reference like this to align all the raw 'reads' to (reads for RNA-seq are often 50-200bp, which is usually unique enough to identify the position this should fit, more on that soon). 
  - GFF3 files, which are not well-standardised and in addition to including start/stop positions for genes and their names, may include all sorts of other metadata.
  - BAM files, which are effectively somewhat standardised and filtered raw reads data. They include positions (start/stop) for ALL the reads, which can be many thousands, and often include a quality score to indicate confidence of alignment and nucleotide along with extra metadata. These files are huge: often gigabytes, and that's just for bacteria...
  - BigWig files can include number of reads per position, which is extracted from the BAM file (in my case, with the 'bamCoverage' program provided by deeptools). A 'bin size' must be provided, which averages across that many values and repeating at every value divisble by that number without a remainder (for example, (1, 2, 3, 4, 5, 6) averaged with a bin size of 3 is (2, 2, 2, 5, 5, 5) output).

JBrowse2 Web was shown to me by post-doc Sam. I used this library to help build the graphs (which reduced my workload, as it is a well-tested genomic browser).
