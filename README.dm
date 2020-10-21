Program use

Overview

BreakAlign is a Perl program downloadable from https://github.com/marchiem/Breakalign along with sample data files necessary
for a basic test run. The program outputs chimaeric sequences aligned to the host genome sequence in what has become a standard format.

Sample command line

> perl breakalign.pl —f chr19_29855576_29856018.fa —r NGSreads.fa —vr LTR.fa

Prior to execution the following key steps are required.

1. BreakAlign will check for a system blastn installation since this dependency is a strict prerequisite for the software to run; if not found,
the user is warned and can provide a directory path where a local version of blastn is downloaded using the -bp switch and write the local blastn
directory path (if the version of blastn predates 2.2.31+ check Options below).

2. FASTA files of the viral LTR sequence and the NGS read file need to be indicated by the —vr and —r switches respectively.

3. Finally, BreakAlign requires either a approximately 500nt human genomic sequence in FASTA format that spans the putative
integration site or a coordinate of such a sequence within the reference human genome sequence. In the sample command line above,
the —f switch directs the program to the former. To use the latter, see Options below.

Options

BreakAlign can retrieve automatically a reference sequence spanning a putative integration site from an installed
copy of the reference human genome. This genome must be downloaded from either the UCSC website (http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/)
or from the NCBI –NIH genomes FTP site (ftp://ftp.ncbi.nih.gov/genomes). The downloads will contain chromosomes as individual FASTA file within a single folder (chr1.fa, chr2.fa, etc). To use this option the following are required.

1. The Bioperl module SeqIO must be installed.

2. Indicate the directory path where the host reference genome sequence is stored using the —fr switch (/path/to/reference_folder/).

3. A coordinate range as a string in ‘chrN:start-end’ format must be provided using the —c switch (e.g. chr16:89577447-89578013).

4. Do not use the —f switch, which overrides this option.
