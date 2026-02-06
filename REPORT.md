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
    

  

Bacteria such as P.aeruginosa are often studied to understand the 

- Learned w
- Sam gave me JBrowse

what ive done:
- w1
	- study rna/dna basic functionality
	- jbrowse is a very 'ideal' visualiser already
		(and just found out its apache2 open src and is ready to be deployed!?)
		
	- day2: more of the same... realised i need to learn a lot of basic defs like RNA transcription from DNA, cDNA reverse transcription from RNA then
	synthesis to produce DNA fragment for rna-seq machine (it only uses DNA, which is more stable).
	Reads are literally just a reading from the rna-sequencer, one run over a particular sample. eg might give: AGCTAGGT. We can use
	this fragment to map it to a gene we know about in the DNA. High read counts might show that the gene is highly expressed.
	
	day3: more in-depth about process. this YT vid gives info. bioinformatics for rnaseq.rnote has my notes from it.
	https://www.youtube.com/watch?v=OxgiyS9Wvww
	
	Phenotypes we might see are result of different RNA levels.
	(phenotypes == total set of observable characteristics, the result of interaction between genetic makeup, aka genotype, and environment factors. eg height, eye colour etc.)
	
	differential expression = how does number of copies of transcript change, and what does it say about the thing.
	EG maybe certain gene more expressed when person sick.
	
	If using JBrowse, likely want:
		- AlignmentsTrack (BAM) → raw reads & splice junctions (ALIGNMENTS - SHOWS WHERE READ IS. MULTIPLE STACKED UP GIVES HIGHER COUNT/EXPRESSION)
		- QuantitativeTrack (BigWig) → the smooth coverage graph (COUNTS SMOOTHED)

=== DAY 4 ===
	day4: need to ask Sam what we really after!
		Options:
		1. I create something from scratch. This has potential to be faster, perhaps better view when zoomed out. But it will be SIMPLE. 
			It's not possible to make something as complex and full-featured as Jbrowse in 9-10 weeks.
		2. We use Jbrowse - and adapt it somehow (but what do we even want?). This may make it harder to change styling or adapt it perfectly
			for what we need but it will have all the features of Jbrowse, and I can make it run offline etc. Should be ways to speed up processing??
			
	Talked with Sam:
		- JBrowse is pretty good. We want it deployable on Daniel's website.
		- Check i can load file types that Sam gave me in ~/Desktop/files_for_jbrowse
		- Main thing is, we want an upload button in addition to standard JBrowse!! This doesn't exist apparently!
		  So we need something where we can upload like a BAM file or whatever (~300mb-400mb!), then the visualisations
		  popup on JBrowse once it loads. This ok - we should be able to add additional stuff to JBrowse so that we can 
		  upload and students in the dept. can all use it.	
		  
	SO - get those files working in JBrowse. 
	I think we'll need to make it so config.json gets changed after an upload, then probs redirects the user to new page with
	sep config.json to show the data that the student uploaded.
	
	
	
	=== DAY 5 ===
	Definitely JBrowse. 
	
	Flow:
	(SERVER)
	Upload data -> check filetypes/structure? -> give appropriate config.json file to render -> redirect onload
	
	Server will need to somehow dynamically serve config.json along w/ appropriate data (via API)
	
	OK - got static HTML doing it ok! VCF and FASTA work. There are some 'gotchas' - see code comments.
	Probably worth talking to Sam again after try out GFF3 files (yikes, those the ones with weird fields).
	
	
WEEK 2
day 1:
	tried to get GFF3 working - seems ok. 
	Realised that we need a bunch of CLI tools (which use a bunch of different things... conda, random installs etc.)
	Decided that using Docker will make things a LOT easier for setup etc.
	
day2:
	got BED working... but it's incredibly slow. Will likely need converted to BigBed first so it's usable.
	BAM is next and last one. still need to figure out how to reconcile 'chromosome1' with 'NC011700.1' with 'NC011700'
	Ok yeah BAM is being a fkn nightmare - sending like 100 requests of 100kb each??? May need j BigWig (essentially BAM w/out reads info).
	/home/eb/Documents/Uni/summer25/pletzerjbrowse/fixed.bam
	Theres just too much I believe. >= 1300 reads on just the one small position i zoomed into? Seems like it should be ok,
	but apparently we might just need to use BigWig here - which is essentially BAM without all the read info (just counts like
	Daniel originally asked for). Time to talk to Sam tomorrow morn...
	
	Points for Sam:
	
	- BAM is too slow (this file at least which has like >1300 reads in small spaces). Is BigWig to show counts good enough if v large? 
	We may potentially be able to LOCALLY serve a FILE instead of using HTTP at all which would likely speed it up, because currently
	its fetching in tiny little chunks of <1Mb... but still...?
	- Another option with BAM: using local files might work better. EG a desktop app that does all this instead of a HTTP server.
	- Most other formats were readable. Some needed conversion first - I can probably try match up mismatched chromosome names
	based on the indexes from the FASTA file.
	
	Got all conversions working (for most part) and saved them in script files / README.md.
	
	Week 3:
	
	Day 1 / 2
	Got search working (jbrowse text-index)
	Got diff colours based on direction of read for GFF genes as sam wanted.
	Customizing the 
	Made checkers for GFF/BED/BAM/FASTA to extract ALL genes. Then we can check for mismatches (finish this later when they know what they want).
	
	Docker setup for deps.
	finished BED -> BigBed. negligible improvement.
	Sam: we only need a few ref FASTA genomes probably - like some combobox to select.
	
	Day 3
	Found out JBrowse has "aliases" file, so, big conversions may be unnecessary and we can just use alias names! THIS would probably be a lot better.
	Apparently performance can be improved via JBrowse Web instead of React component as it utilizes workers etc. better.
	May be more annoying to customise potentially.
	SO FAR - quite a LOT FASTER with full 'JBrowse Web'!! LIKE WAY BETTER.
	Did some theming!
	OK - converted to JBrowse Web version. Its quite obviously faster - doesnt freeze when scrolling. Will keep with this version and have to conv
	some stuff over to it. MUCH better (uses web workers so it doesnt freeze the main React thread)
	
	
Week 4:

Sam wants:
	- see if we can get diff colours for diff strains or something.
	- script to warn against mismatch (mostly done)

Daniel wants:
	- test my JBrowse site with the data from that other website (make sure they look similar)
	- Daniel's concern about the BigWig file looking different to the other website (I BELIEVE ITS DUE TO MY DATA PROCESSING. if you zoom in - its
	averaged over 50 indices! this is because in bamCoverage, the default binSize is 50. we likely want lower then!!)
	- try log10 transform
	- no uploads, just change between datasets.
	- this deployed to wordpress site (static site will make this easy i think).
