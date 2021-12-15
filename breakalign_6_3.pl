#!/usr/local/bin/perl -w
use Data::Dumper;
use strict;
use Cwd;
use File::Copy;
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path);

my $br_path = abs_path(__FILE__); #this is to get the breakalign path, needed to call it within the code
#die $br_path; it would show this /Users/emarchi/Dropbox/scripts/breakalign_6_1.pl

my $help = 0;
my $uref = "";
my $coo = ""; # coordinates can be given in format chrN:start-end
my $hg19_ref = ""; # can be given; will probably be in form /databank/raw/hg19_full/
my $vir_ref = ""; # can put in path to LTR file here rather than in command line
my $reads = ""; # reads in fasta format in input
my $word_s = '10';
my $word_s2 = '20';
my $blast_path = ""; # can put in path to BLAST here rather than in command line; will probably be like /.../ncbi-blast-2.2.29+/bin/
#my $k_f = ""; #option to keep temporary files
my $k_f = 0; #option to keep temporary files, it seems to work better initializing the variable with 0 rather than empty character
my $bed_co = ""; #file with multiple coordinates in bed format to test 
my $k_d = 0; #option to keep just the database of all reads in input, to avoid to build it for each site tested
my $n_t = 1; #number of threads
#die $help;
unless ($blast_path) {
	$blast_path = qx(which blastn);
	$blast_path =~  s/blastn$//;
	chomp($blast_path); # remove new line character if we get blast path from 'which' system command
	$blast_path =~ s/\r\n$//; # just in case chomp does not work
	$blast_path =~ s/\r$//;# just in case...

}

unless ($blast_path) { die "\nERROR: blastn not found"; }
else { print "NCBI Blast installation detected\n"}

system($blast_path."blastn -version ") == 0 or die "I cannot find blastn at $blast_path";

GetOptions ("help"  => \$help, "f=s" => \$uref, "c=s" => \$coo, "fr=s" => \$hg19_ref, "r=s" => \$reads,
            "ws=s" =>\$word_s, "ws2=s" => \$word_s2, "bp=s" => \$blast_path, "vr=s" => \$vir_ref, "kt" => \$k_f,"bf=s" => \$bed_co,"kd" => \$k_d,
						"nt=i" =>\$n_t)
or die("Error in command line arguments\n");

if ($help)
{
	  
    print "usage: perl path/to/Breakalign.pl [options]\n";
		print "ESSENTIAL\n";
		print "-r   FILE    File containing sequences reads\n";
		print "-vr  FILE    File containing LTR sequence reference\n";
		print "AND EITHER\n";
		print "-f   FILE    ~500nts of human reference sequence spanning integration\n";
		print "OR\n";
		print "-fr  PATH  Reference folder containing FASTA files for each chromosome\n";
		print "     AND EITHER\n";
		print "     -c   STR   Genomic coordinates <Chr>:<Start>-<Stop> (e.g. chr16:89577447-89578013)\n";
		print "     OR\n";
		print "     -bf  FILE  Multiple coordinates in bed format\n\n";
		print "OPTIONAL\n";
    print "-h		      This menu\n";
		print "-bp  PATH  Path to alternative blastn installation\n";
		print "-kt        Keep all temporary files (to be used for a single site inspection)\n";
		print "-kd        Keep reads database files (to be used for large input reads file)\n";
		print "-nt  INT   Number of threads to use in blastn query command\n";
    print "-ws  INT   Minumun word size alignment of reads to the reference sequence [10]\n";
    print "-ws2 INT   Minumun word size alignment of reads to the LTR sequence [10]\n";
    
    
		
		
		
		#die $help; supposed to be 1 when help option is chose. ok
    exit;
}

my $help_out = "Breakalign requires the following switches:\n";
$help_out .= "-r          NGS reads file\n";
$help_out .= "-vr         LTR sequence file\n";
$help_out .= "AND EITHER\n";
$help_out .= "-f          File containing ~500nts of human reference sequence spanning integration\n";
$help_out .= "OR\n";
$help_out .= "-fr         Path to reference folder containing fasta files for all chromosomes\n";
$help_out .= "            AND EITHER\n";
$help_out .= "            -c         Genomic coordinate in form <Chr>:<Start>-<Stop>\n";
$help_out .= "            OR\n";
$help_out .= "            -bf        Multiple coordinates in bed file\n\n";
$help_out .= "Check usage with Breakalign.pl -h\n";
if ($vir_ref ne "") {print "LTR reference sequence in folder: ".$vir_ref."\n"} else {print "Please provide LTR reference sequence\n\n".$help_out;exit;};
if ($reads ne "") {print "Reads to align to reference: ".$reads."\n"} else {print "Please provide reads to align to the references\n\n".$help_out;exit;};
if ($uref ne "")
{
    print "I will use user provided fasta file as reference: $uref.\n"
    
}
elsif ($coo)
    {
        if ($uref eq "") {require Bio::SeqIO};    
        if ($coo ne "") {print "Reference sequence will be extracted from coordinates ".$coo."\n"} else {print "Please provide Genomic coordinates\n\n".$help_out;exit;};
        if ($hg19_ref ne "") {print "from raw fasta reference in folder: ".$hg19_ref."\n"} else {print "Please provide path to folder with reference fasta files\n\n".$help_out;exit};
        #print "LTR reference sequence given: ".$vir_ref."\n";#just printed above and not meant to be here
        #print "Reads to align to reference given: ".$reads."\n";#just printed above
    }
elsif ($bed_co)
		{
			print "Coordinates provided in bed file format\n";
		}
else
{
  print "Please provide reference file and suitable input files\n\n".$help_out;exit;  
}

# dealing with candidate region genome coordinates, input must be like chr3:76076289-76076636

my $fasta_ref;

unless ($uref)
{
		if ($coo) #in case coordinates are provided
		{
			my @ar_c = split(":",$coo);
			my $chr_c = $ar_c[0];
			my @st_en = split("-",$ar_c[1]);
			my $st = $st_en[0];
			my $en = $st_en[1];


			$fasta_ref = "hg19_".$chr_c."_".$st."_".$en.".fa";
			open(FASTAREF,">$fasta_ref") or die "I cannot open $fasta_ref";
			print FASTAREF ">hg19_".$chr_c."_".$st."_".$en."\n";
			$hg19_ref.="/" if(not $hg19_ref =~/\/$/ ); #this will fix the missing slash in the path
			my $seqio = Bio::SeqIO->new(-file => $hg19_ref.$chr_c.".fa", -format => "fasta");

			while(my $seq = $seqio->next_seq)
			{ 
			    print FASTAREF uc($seq->subseq($st,$en)), "\n";
			    print uc($seq->subseq($st,$en)), "\n";
			}

			close (FASTAREF);
		}
		elsif($bed_co) #in this case we use a bed file with set of coordinates
		{
			#I call breakalign inside breakalign, in loop, with input coordinatesfrom each line of the bed file
			open(BCO,"$bed_co") or die "I cannot open "."$bed_co";
			while(<BCO>)
			{
				my $line_co = $_;
			  chomp($line_co);
			  unless($line_co) {next;}

				my @co = split /\s+/, $line_co;
				
				unless ($co[0] =~ m/^chr|^Chr/) {$co[0] = "chr".$co[0]};
				if ($co[0] =~ m/^Chr/) {$co[0] =~ s/^C/c/};
				$coo = $co[0].":".$co[1]."-".$co[2]; #Recreating the correct input format for breakalign
				
				# input and output folder
				# I am reapeting this later maybe it should be done once at the beginning of the code
				my $input_dir = getcwd;# I just need this
				my $output_dir = getcwd;#here I do not need this? 

				#if ($input_dir =~ m/ /) {$input_dir =~ s/ /\\ /};# I was doing this for windows system?
				#if ($output_dir =~ m/ /) {$output_dir =~ s/ /\\ /};

				my $input_dir_b = $input_dir; #I am not sure why we had to do the following, it was for a bug

				if ($input_dir =~ m/\\ /) {$input_dir_b =~ s/\\ / /;};
				if (-e $input_dir_b."/".$bed_co)
				{
					# make the database using the input reads, and it will be done only once
					
					my $cmd_br = "perl ".$br_path." -fr ".$hg19_ref." -vr ".$vir_ref." -c ".$coo." -r ".$input_dir."/".$reads." -ws ".$word_s." -ws2 ".$word_s2." -kd"." -nt ".$n_t;
					#print "executing: $cmd_2\n";
					#print FILER "executing: $cmd_br\n";
					system($cmd_br) == 0 or die "system cmd [$cmd_br] failed ($?): $!";
            
         }
				
			}
			exit; #once we ran breakling for each line coordinate we can exit
			
		}
}
else
{
    $fasta_ref = basename($uref);
}

# dealing with fasta file of read sequences; possibly not all of them, only the ones mapping the candidate region or known to contain the LTR

my $reads_f = basename($reads);
my $word_size = $word_s;
my $word_size_b = $word_s2;

# labelling for unique output

my $label;
unless ($uref)
{
    $label = $fasta_ref;
    $label =~  s/\.[^.]+$//;
}
else {$label ="User_reference"}

# generating a new input file with unique and numeric ids

my $new_reads_f= "Unique_".$reads_f;

open(FH, '>', $new_reads_f) or die "Could not open file '$new_reads_f' $!";

open(FIN, $reads) or die "Could not open file '$reads' $!";

my $n_line = 0;
while (<FIN>)
{
    my $line = $_;
    chomp($line);
    unless($line) {next;}
    
    $line =~ s/\r\n$//; # just in case chomp would not work
    $line =~ s/\r$//;# just in case... to deal with any os
    if ($line =~ /^>/)
    {
        $n_line++;
        $line = ">".$n_line;
    }
    print FH $line."\n";
}


$reads_f = $new_reads_f;

# input and output folder

my $input_dir = getcwd;
my $output_dir = getcwd;

if ($input_dir =~ m/ /) {$input_dir =~ s/ /\\ /};
if ($output_dir =~ m/ /) {$output_dir =~ s/ /\\ /};

# copy the reference sequence in the input dir in case it is not already in the working directory

#copy($uref,$input_dir);

# to handle escape character in -e or open in perl

my $output_dir_b =  $output_dir;
if ($output_dir_b  =~ m/\\ /) {$output_dir_b =~ s/\\ / /;}

# my temp output

my $temp_out = $output_dir."/temp_br_al";

# my report output

my $rep_f = $output_dir."/Report_".time().".txt"; # label temporary the main output with a time stamp

my $out_file = $output_dir."/align_for_breakpoint"."_".$label.".txt";

# make blast output file

my $blast_f = $temp_out."/blast_out_ref_";

# breakalign fasta input

my $br_input =$temp_out."/hits_seq_".$label.".fa";
my $br_input_2 =$temp_out."/hits_seq_K03455_".$label.".fa";

# open output files

if ($rep_f =~ m/\\ /) {$rep_f =~ s/\\ / /};
open FILER, ">$rep_f" or die "unable to open $rep_f $!";

if ($out_file =~ m/\\ /) {$out_file =~ s/\\ / /};
open(OUT_F,">$out_file") or die "I cannot open $out_file\r\n" ;

# create a temp folder if there is not one

if ($temp_out =~ m/\\ /) {$temp_out =~ s/\\ / /}; # because -e does not want escape character
unless (-e $temp_out)
{
    my $cmd_1 = "/bin/mkdir $temp_out";
    #print "executing: $cmd_1\n";
    print FILER "executing: $cmd_1\n";
    system($cmd_1) == 0 or die "system cmd [$cmd_1] failed ($?): $!";  # fork, execute $cmd and return exit status
}

# running makeblastdb to make the blast database of the candidate region in the correct format 

my $input_dir_b = $input_dir;

if ($input_dir =~ m/\\ /) {$input_dir_b =~ s/\\ / /;};
unless (-e $input_dir_b."/blastdb_reads.nin" or -e $input_dir_b."/blastdb_reads.00.nin" or -e $input_dir_b."/All_reads_db.00.nhr" or -e $input_dir_b."/All_reads_db.nin")
{

    if (-e $input_dir_b."/".$reads_f)
    {
        # make the database using the reads, not the genomic region
        
        my $cmd_2 = $blast_path."/makeblastdb -in ".$input_dir."/".$reads_f." -dbtype nucl -parse_seqids -out ".$output_dir."/All_reads_db -title 'All_reads_db'";
        #print "executing: $cmd_2\n";
        print FILER "executing: $cmd_2\n";
        system($cmd_2) == 0 or die "system cmd [$cmd_2] failed ($?): $!";
            
    }
    else
    {
        print "The reference file $reads_f does not exist or it is not in the input directory: ".$input_dir;
    }
}
else {print "Blast database blastdb_reads exists in the current directory, it will not be redone\n"}

# running blast using the reads sequences in input versus the candidate region

my $blast_f_label =  $blast_f.$label;
if ( $blast_f_label =~ m/\\ /) {$blast_f_label =~ s/\\ / /};

unless (-e $blast_f_label)
{
    my $cmd_3 = $blast_path."/blastn -query ".$input_dir."/".$fasta_ref." -num_threads ".$n_t." -db 'All_reads_db' -word_size $word_size -outfmt '7 sseqid qstart qend sstart send qseq sseq evalue bitscore score length qframe sframe' -out ".$blast_f_label; 
    #print "executing: $cmd_3\n";
    print FILER "executing: $cmd_3\n";
    system($cmd_3) == 0 or die "system cmd [$cmd_3] failed ($?): $!";

}
else
{
    print "The genomic region in file $label has been analyzed with blastn before."."\n"."The output file is in ".$blast_f_label."\n"; 
}

####################################################################################
# converting the reference into an array and then into a hash; able to handle indels
####################################################################################

open(FA, $fasta_ref) or die "I cannot open $fasta_ref";
my @array_ref;
my %hash_ref;
my $ref;
while(<FA>)
{
    my $line = $_;
    chomp($line);
    $line =~ s/\r\n$//; # just in case chomp would not work
    $line =~ s/\r$//;# just in case...
    if ($line =~ /^>/) {next};
    $ref .= $line; # it should handle new lines characters in the reference sequence 
}

close(FA);

@array_ref = split('',$ref);

my $i =1 ;
foreach (@array_ref)
{
    $hash_ref{$i}{$_} = ''; # so each position can be either a base or a gap '-'
    $i++;
    
};
####################################################################################
# storing reads with hits into a hash; hits in the region or the fasta file provided
####################################################################################

# not scanning for the significant read sequences in the file in input as it could be too big; use blastdbcmd instead to retrieve the sequence with hits
# checking what read ids have significant hits and keep the sequence

my %hits;
my %query_al;
my %read_al;

if ($blast_f_label =~ m/\\ /) {$blast_f_label =~ s/\\ / /};
my %seen;
open(BLASTF,$blast_f_label) or die "I cannot open $blast_f_label\n";
while (<BLASTF>)
{
    if ($_ =~ /^#/) {next;}
    my $line = $_;
    chomp($line);
    my @ar = split('\t',$line);
    # get only the best hit of each alignment
    if (exists($seen{$ar[0]})){next;}
    
    if ($ar[3] < $ar[4])
    {
        $hits{$ar[0]}{'sense'} = 'F'; # the alignment is forward
        $hits{$ar[0]}{'st_al'} = $ar[1];
        $hits{$ar[0]}{'st_al_read'} = $ar[3];
        
    }
    else
    {
        $hits{$ar[0]}{'sense'} = 'R'; # the alignment is with the reverse complement of the read
        $hits{$ar[0]}{'st_al'} = $ar[1];
        $hits{$ar[0]}{'end_al_read'} = $ar[3];
       
    }
    my $align_q = $ar[5];
    my @ar2 = split('',$align_q);
    my $p = $ar[1]; # alignment on the query sequence
    foreach (@ar2)
    {
           my $base = $_;
           $query_al{$ar[0]}{$p}{$base} = '';
           $p++;
    }
    
    my $align_r = $ar[6];
    my @ar3 = split('',$align_r);
    my $p2 = $ar[1]; # put the read alignment in the query position coordinates
    foreach (@ar3)

    {
       my $base = $_;
       $read_al{$ar[0]}{$p2}{$base} = '';
       $p2++;
    }
    
    # now get the sequence for each read id using blastdbcmd

    my $cmd_4 = $blast_path."/blastdbcmd -entry 'lcl|".$ar[0]."' -db 'All_reads_db' >> ".$br_input;
    #print "executing: $cmd_4\n";
    print FILER "executing: $cmd_4\n";
    system($cmd_4) == 0 or die "system cmd [$cmd_4] failed ($?): $!";
    $seen{$ar[0]} = '';
    
}
close(BLASTF);

#1. Get the sequences with significant hits
#2. blast the sequences with hits versus K03455 (making a db of them)
#3. Producing a second filtered input for breakalign with only reads having hits with K03455, likely to contain the LTR
# Important: we have to make sure that each sequence has a unique id. Major difference with previous version in openning the fasta file with sequences retrived  by blastdbcmd

my $br_input_b =  $br_input;
if ($br_input_b  =~ m/\\ /) {$br_input_b =~ s/\\ / /;}
open(RE,"$br_input_b") or die "I cannot open "."$br_input_b";

# first remove the new lines; the script below expects fasta format with no new line characters

my $conc_line;
my $firstline = 1;
my $outvar;
my $IDD;
while(<RE>)
{
    my $line = $_;
    chomp($line);
    $_ =~ s/\n$//; # in case chomp would not work...
    $line =~ s/\r$//; # remember to remove \r too in case on a mav
    $line =~ s/\r\n$//; # not needed but one or the other should work
    if ($line =~ /^>/)
    {
        if ($firstline == 0)
        {
            $outvar.= $IDD."\n".$conc_line."\n";
        }
        $IDD = $line;
        $conc_line = '';
        $firstline = 0;
        next;
    };
    $conc_line .= $line; # should handle new lines characters in the reference sequence 
    
}
$outvar.= $IDD."\n".$conc_line."\n"; # this is to print the last line
close(RE);
open(REWR, ">".$br_input_b) or die "I cannot write on NEW_".$br_input_b."\n";
print REWR $outvar;
close(REWR);

####################################################################################
# finding the hits for the LTR in the mapping reads
####################################################################################

if (-e $br_input_b)
    {
        # make the database using the reads, not the genomic region
        my $cmd_5 = $blast_path."/makeblastdb -in ".$br_input_b." -dbtype nucl -parse_seqids -out ./hits_reads_db_".$label." -title 'hits_reads_db_".$label."'";
        #print "executing: $cmd_5\n";
        print FILER "executing: $cmd_5\n";
        system($cmd_5) == 0 or die "system cmd [$cmd_5] failed ($?): $!";
            
    }
    else
    {
        print "The fasta file $br_input with read sequences having hits to the candidate regions does not exist";
    }

if (-e $output_dir_b."/hits_reads_db_".$label.".nhr")
{
    my $cmd_6 = $blast_path."/blastn -query ".$vir_ref." -num_threads ".$n_t." -db 'hits_reads_db_".$label."' -word_size ".$word_size_b." -outfmt '7 sseqid qstart qend sstart send qseq sseq evalue bitscore score length qframe sframe' -out ".$temp_out."/K03455_blast_out_".$label;
    #print "executing: $cmd_6\n";
    print FILER "executing: $cmd_6\n";
    system($cmd_6) == 0 or die "system cmd [$cmd_6] failed ($?): $!";

}
else
{
    die "The second database of reads mapping the genomic region ".$output_dir_b."/hits_reads_db_".$label." file does not exist or it is not in the correct directory: ".$output_dir_b;
}
# generating the second breakalign input with only reads having hits for K03455

my $temp_out_b =  $temp_out;
if ($temp_out_b  =~ m/\\ /) {$temp_out_b =~ s/\\ / /;}


open(BLASTH,$temp_out_b."/K03455_blast_out_".$label) or die "I cannot open $temp_out_b"."/K03455_blast_out_".$label;
while (<BLASTH>)
{
    if ($_ =~ /^#/) {next;}
    my $line = $_;
    chomp($line);
    my @ar = split('\t',$line);
    
    my $cmd_7 = $blast_path."/blastdbcmd -entry 'lcl|".$ar[0]."' -db 'hits_reads_db_".$label."' >> ".$br_input_2;
    #print "executing: $cmd_7\n";
    print FILER "executing: $cmd_7\n";
    system($cmd_7) == 0 or die "system cmd [$cmd_7] failed ($?): $!";
}
close(BLASTH);

# first remove the new lines; the script below expects fasta format with no new line characters

my $br_input_2_b =  $br_input_2;
if ($br_input_2_b  =~ m/\\ /) {$br_input_2_b =~ s/\\ / /;}


if (-e $br_input_2_b)
{
    open(RE2,"$br_input_2_b") or die "I cannot open "."$br_input_2_b";
    my $conc_line_2;
    my $firstline_2 = 1;
    my $outvar_2;
    my $IDD_2;
    while(<RE2>)
    {
        my $line = $_;
        chomp($line);
        $_ =~ s/\n$//; # in case chomp would not work...
        $line =~ s/\r$//; # remember to remove \r too in case we are in a mac
        $line =~ s/\r\n$//; # not needed but one or the other should work
        if ($line =~ /^>/)
        {
            if ($firstline_2 == 0)
            {
                $outvar_2.= $IDD_2."\n".$conc_line_2."\n";
            }
            $IDD_2 = $line;
            $conc_line_2 = '';
            $firstline_2 = 0;
            next;
        };
        $conc_line_2 .= $line; # should handle new line characters in the reference sequence 
        
    }
    $outvar_2.= $IDD_2."\n".$conc_line_2."\n"; # this is to print the last line
    close(RE2);
    open(REWR2, ">"."$br_input_2_b") or die "I cannot write on NEW "."$br_input_2_b"."\n";
    print REWR2 $outvar_2;
    close(REWR2);
# end of second input to breakalign
}
else
{
	print "File: ".$br_input_2_b." is missing or not gerated,\n so the reads do not aligned to virus/retrotrasposon reference\n";
	print "Temporary files will be deleted"."\n";
	Deletefile();
	
	exit;
}


##########################################################################
# get read ids and sequences and make the reverse complement sequences; keep only the reads having hits to the virus sequence
##########################################################################

open(RE,$br_input_2_b) or die "I cannot open $br_input_2_b\n"; # $br_input_2_b contains reads found to have hits to the region and the virus at the same time
my $r_id;
my %h_reads_s; # hashes with the read sequences saved position by position, like the reference
my %h_reads_rc; # not used at present but might help solve the problem with bugs for too long printing
    
my $length_re; # important to work out the proper alignment
while(<RE>)
{
    my $line = $_;
    chomp($line);
    $_ =~ s/\n$//; # in case chomp would not work
    $line =~ s/\r$//; #remember to remove \r too in case we are in a mac
    $line =~ s/\r\n$//; # not needed but one or the other should work
    if ($line =~ /^>/)
        {
            if ($line =~ /(\d+)\s/)
            {
            $r_id = $1;
            }
        }
        else
        {
            my $seq = $line;
            if (exists ($hits{$r_id})) # store only read sequences having significant hits
            {
                $length_re = 0; # get the length of the read sequence
                my @ar = split('',$seq);
                foreach (@ar)
                {
                    $length_re++;     
                }
                            
                my $rsir;
                if ($hits{$r_id}{'sense'} eq 'F')
                {
                    my $rsir = ($hits{$r_id}{'st_al'} - $hits{$r_id}{'st_al_read'})  + 1; # read start in reference
                    my $i=$rsir;
                    $length_re = 1; # work out the length of the read
                    foreach (@ar)
                    {
                        $h_reads_s{$r_id}{$seq}{$i}{$_}= $_;
                        $i++;             
                    }
                                    
                }
                if ($hits{$r_id}{'sense'} eq 'R')
                {
                    my $revc = reverse_complement($seq);
                    my @ar_rc = split('',$revc);
                    my $rsir = $hits{$r_id}{'st_al'} - ($length_re - $hits{$r_id}{'end_al_read'});
                    my $ii=$rsir;
                    foreach (@ar_rc)
                    {
                        $h_reads_s{$r_id}{$revc}{$ii}{$_}= $_;
                        $ii++;
                     }
                }
            }   
        }
    }

close(RE);

#########################################################
# print out the alignment 
#########################################################

# the genomic region first
# make the header with the reference sequence
foreach my $pos ( sort _numerically_asc keys %hash_ref)
{
    foreach my $b (keys %{$hash_ref{$pos}})
    {
        print OUT_F $b;
    }
}

print OUT_F "\n"; 

my $out; # output variable to print out
my $indent;
foreach my $ids (sort keys %hits)
{
    unless (exists($h_reads_s{$ids})){next;}  # print out the reads found to have hits

    foreach my $seq (keys %{$h_reads_s{$ids}}) # get the indentation using the query alignment start, which is the start of the alignment on the genomic region
    {
        foreach my $pos (sort _numerically_asc keys %{$h_reads_s{$ids}{$seq}}) 
        {
            $indent = $pos;
        last;
        }
    }
    
    my $line_out .= ' ' x ($indent - 1 );
    my $inser = '';
    my $dele='';  
    foreach my $seq (keys %{$h_reads_s{$ids}}) # check which one of the four possible sequences for each id to print (the mate that aligns and in which direction)
    {
         foreach my $pos (sort _numerically_asc keys %{$h_reads_s{$ids}{$seq}}) 
        {
            unless (exists ($query_al{$ids}{$pos}))
            {
                foreach my $b (keys %{$h_reads_s{$ids}{$seq}{$pos}}) # take the only base possible there
                {
                    $line_out .= lc($b);
                }              
            }
            else
            {
                foreach my $b (keys %{$query_al{$ids}{$pos}}) # takes the only base possible there
                {
                    #unless($b eq '-'){$out .= uc($b)};
                    if($b eq '-')
                    {
                      $line_out =~ s/^\s//;# remove a white space from the indentation
                      $inser.=(keys %{$read_al{$ids}{$pos}})[0]; # in the same position I check what the base should be
                      $b = (keys %{$read_al{$ids}{$pos}})[0]; # prints the mismatched base in read alignment
                      $line_out .= lc($b); # mark and make visible the insertion with a lower case letter
                      next; # this is just to the next base; otherwise a line below would make an upper case again
                    } 
                    if((keys %{$read_al{$ids}{$pos}})[0] eq '-') # checks if there is deletion -- would be showed in the read alignment with a "-"
                    {
                      $dele.=(keys %{$query_al{$ids}{$pos}})[0]; # in the same position check what base should be
                      
                      $b = (keys %{$read_al{$ids}{$pos}})[0] # same as above change, but I would print the "-" instead of the base (from query alignent anyway) because in the read it is missing anyway
                    }
                    $line_out .= uc($b); # the "-" character is not affected
                }        
            }
        }
    }
    $out .= $line_out;
    $out .= "   ".$ids;
    if ($inser ne ''){ $out .= " INSERTION: ".lc($inser)." ";}
    if ($dele ne ''){$out .= " DELETION: ".$dele." ";}
    $inser = '';
    $dele = '';
    $out .="\n";
       
}

print OUT_F $out;
close (OUT_F);

########################
# delete temporary files
########################
#Usinga subroutine now
 
unless ($k_f){Deletefile();}

#unless ($k_f)
#{
#  # remove temp files (in version 5.6 the temp directory with hits was removed)
#  print "Deleted temporary files\n";
#  opendir(DIR, $temp_out);
#  my @t_files = readdir DIR;
#  
#  foreach my $t_file (@t_files)
#  {
#   unlink($temp_out."/".$t_file);
#  #print "Deleted temporary file: $t_file\n";
#  }
#  closedir(DIR);
#  rmdir($temp_out)  or warn "couldn't remove $temp_out: $!";
#
#  #opendir(DIR, $input_dir);
#  #my @files = grep(/^All_reads_|^hits_reads_/,readdir(DIR));
#  #foreach my $file (@files)
#  #{
#  # unlink($input_dir."/".$file);
#  # #print "Deleted reads database file: $file\n";
#  #}
#  #closedir(DIR);
#	
#	opendir(DIR, $input_dir);
#  my @files_h = grep(/^hits_reads_/,readdir(DIR));
#  foreach my $file (@files_h)
#  {
#   unlink($input_dir."/".$file);
#   #print "Deleted database file of reads aligning to Virus: $file\n";
#  }
#  closedir(DIR);
#	
#	unless ($k_d)
#	{
#		opendir(DIR, $input_dir);
#		my @files_A = grep(/^All_reads_/,readdir(DIR));
#		foreach my $file (@files_A)
#		{
#		 unlink($input_dir."/".$file);
#		 #print "Deleted reads database file: $file\n";
#		}
#		closedir(DIR);
#	}
#}

##############################################################
# make the html output by reusing br_text_to_html_mark.pl

br_txt_to_html($out_file);

#####################################
# subroutines                     ###
#####################################

sub reverse_complement {

        my $dna = shift;

	# reverse the DNA sequence

        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence

        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}



sub _numerically_asc

{
	$a <=> $b;
}



sub _numerically_dis

{
	$b <=> $a;
}

sub br_txt_to_html
{
     my $br_out_txt = shift;
     my $br_out_html = $br_out_txt;
     $br_out_html =~ s/txt$/html/;
     
     my $blast_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch&USER_FORMAT_DEFAULTS=on&SET_SAVED_SEARCH=true&PAGE=MegaBlast&PROGRAM=blastn&QUERY=";
     
     open(BR_TXT, $br_out_txt) or die "I cannot open $br_out_txt";
     open(BR_HTML, ">$br_out_html") or die "I cannot open $br_out_html $!";
     
     my $indent ='';
     my $r_id;
     my $html_out="<html>
                     <head>
                     <style>
                     mark { 
                           background-color: white;
                           color: red;
                          }
                     </style>
                    </head>
                    <p style=\"display:inline; line-height:0.8; white-space: nowrap; \"> <font face=\"Courier New\" size=\"1\">";
    while (my $line = <BR_TXT>)
    {
        chomp $line;
        
        if ($line =~ /^[A-Z]/) # must be only the first line, the reference
        {
            $html_out.= $line."<br>";
            next;
        }
        
        if ($line =~ /(^\s+)/)
        {
           
            $indent = $1;
            $line =~ s/^\s+//;
            
            if ($line =~ /(\s+\d+)/) # the id is the only numeric and not the last
            {
                $r_id = $1;
                $line =~ s/\s+\d+//;  
                if ($line =~ m/(\sINSERTION:\s[a-zA-Z]+|\sDELETION:\s[a-zA-Z]+)/)
                {
                  $r_id .= $1;
                  $line =~ s/\sINSERTION:\s[a-zA-Z]+\s|DELETION:\s[a-zA-Z]+//;
                }
            }
            
            if ($line =~ /(^[\p{Lowercase}]{8,})/)
            {
                my $vr_l = $1;
                $line =~ s/^[a-z]{8,}//;
                $html_out.="&nbsp;" x length($indent) ."<mark><a href="."$blast_url"."$vr_l"." style=\"color: Red\" >"."$vr_l"."</a></mark>".$line.$r_id."<br>\n";
            }
            elsif($line =~ /([a-z]{8,})/)
            {
            
                my $vr_r = $1;
                $line =~ s/[a-z]{8,}//;
                $html_out.="&nbsp;" x length($indent) .$line."<mark><a href="."$blast_url"."$vr_r"." style=\"color: Red\" >"."$vr_r"."</a></mark>".$r_id."<br>\n";
            }
        }
    }
    $html_out.= "</font></p></html>";
    print BR_HTML $html_out;
    close(BR_TXT);
    close(BR_HTML);
        
}
##########################
#Delete files suboroutine#
##########################
sub Deletefile
{
		print "Deleted temporary files\n";
		opendir(DIR, $temp_out);
		my @t_files = readdir DIR;
		
		foreach my $t_file (@t_files)
		{
		 unlink($temp_out."/".$t_file);
		#print "Deleted temporary file: $t_file\n";
		}
		closedir(DIR);
		rmdir($temp_out)  or warn "couldn't remove $temp_out: $!";
	
		#opendir(DIR, $input_dir);
		#my @files = grep(/^All_reads_|^hits_reads_/,readdir(DIR));
		#foreach my $file (@files)
		#{
		# unlink($input_dir."/".$file);
		# #print "Deleted reads database file: $file\n";
		#}
		#closedir(DIR);
		
		opendir(DIR, $input_dir);
		my @files_h = grep(/^hits_reads_/,readdir(DIR));
		foreach my $file (@files_h)
		{
		 unlink($input_dir."/".$file);
		 #print "Deleted database file of reads aligning to Virus: $file\n";
		}
		closedir(DIR);
		
		unless ($k_d)
		{
			opendir(DIR, $input_dir);
			my @files_A = grep(/^All_reads_/,readdir(DIR));
			foreach my $file (@files_A)
			{
			 unlink($input_dir."/".$file);
			 #print "Deleted reads database file: $file\n";
			}
			closedir(DIR);
		}
}