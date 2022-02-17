Program use

Overview

BreakAlign is a Perl program downloadable from https://github.com/marchiem/Breakalign along with sample data files necessary for a simple test run. The program outputs chimaeric sequences aligned to the host genome sequence in what has become a standard format. Typing breakalign.pl -help will show the options.

Sample command lines (note, always use hyphen not dash symbols and ensure no whitespaces in directory names)

> perl breakalign.pl -f chr19_29855576_29856018.fa -r NGSreads.fa -vr LTR.fa

> perl breakalign.pl -f chr21_14511311_14511999.fa -r NGSreads.fa -vr LTR-HIV.fa


Prior to execution the following key steps are required.

1. BreakAlign will check for a system blastn installation since this dependency is a strict prerequisite for the software to run; if not found, the user is warned and can provide a directory path where a local version of blastn is downloaded using either the -bp switch or inserting into the script at line 20 (if the version of blastn predates 2.2.31+ it should still run but we cannot guarantee that).

2. FASTA files of the viral LTR sequence and the NGS read file need to be indicated by the -vr and -r switches respectively.

3. Finally, BreakAlign requires either an approximately 500nt human genomic sequence in FASTA format that spans the putative integration site or a coordinate of such a sequence within the reference human genome sequence. In the sample command line above, the -f switch directs the program to the former. To use the latter, see Options below.

Options

BreakAlign can retrieve automatically a reference sequence spanning a putative integration site from an installed copy of the reference human genome. This genome must be downloaded from either the UCSC website (http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/) or from the NCBI–NIH genomes FTP site (ftp://ftp.ncbi.nih.gov/genomes). The downloads will contain chromosomes as individual FASTA file within a single folder (chr1.fa, chr2.fa, etc). To use this option the following are required.

1. The Bioperl module SeqIO must be installed.

2. Indicate the directory path where the host reference genome sequence is stored using the -fr switch (/path/to/reference_folder/).

3. And either
	3a. a coordinate range as a string in ‘chrN:start-end’ format can be provided using the -c switch (e.g. chr16:89577447-89578013).
	3b.	a bed file with multiple coordinates can be provided using the -bf switch (e.g. coord.bed). Sample command with toy chromosome sequences provided is > perl breakalign.pl -bf coord.bed -r NGSreads.fa -vr LTR.fa -fr Human_ref/

4. Do not use the -f switch, which overrides this option.

5. The -kd switch will keep and reuse the initial blastn database, which is the most time-consuming part of BreakAlign (remember to delete this manually after you have finished or it will be reused).