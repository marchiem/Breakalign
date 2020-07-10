#!/usr/local/bin/perl -w
use Data::Dumper;
use strict;
use Cwd;
use Bio::SeqIO;
use File::Copy;

use Getopt::Long;
use File::Basename;

my $help = "";
my $uref = "";
my $coo = ""; #coordinates to be given in format chrN:start-end
my $hg19_ref = "/databank/raw/hg19_full/"; # default value
my $vir_ref = "/home/manu/HIV_ref/K03455_LEFT_LTR.fasta";
my $reads = ""; # reads in fasta format in input
my $word_s = '10';
my $word_s2 = '20';
my $blast_path = qx(which blastn);
my $k_f;

#die $blast_path;

$blast_path =~  s/blastn$//;
#die $blast_path; #/usr/bin/

chomp($blast_path); #remove new line character if we get blast path from 'which' system command
$blast_path =~ s/\r\n$//; #just in case chomp would not work
$blast_path =~ s/\r$//;#just in case...


GetOptions ("help"  => \$help, "f=s" => \$uref, "c=s" => \$coo, "fr=s" => \$hg19_ref, "r=s" => \$reads,
            "ws=s" =>\$word_s, "ws2=s" => \$word_s2, "bp=s" => \$blast_path, "vr=s" => \$vir_ref, "kt" => \$k_f )   # 
or die("Error in command line arguments\n");

unless ($blast_path) { die "\nERROR: blastn not found"; }
else { print "NCBI Blast installation detected\n"}

system($blast_path."blastn -version ") == 0 or die "I cannot find blastn at $blast_path";


if ($help)
{
    print "usage: perl path/to/Breakalign.pl [options] -r <in.fasta>\n";
    print "-help        This menu\n";
    print "-f   FILE    User provided reference suquence\n";
    #print "-f   FILE    User provided reference suquence\n";
    print "-fr          Path to reference folder cointaining fasta files for each chromosome\n";
    print "-vr          path to LTR sequence reference\n";
    print "-c   STR     Genomic coordinates <Chr>:<Start>-<Stop> (e.g. chr16:89577447-89578013)\n";
    print "-r   FILE    File containing sequences reads\n";
    print "-ws  INT     minumun word size alignment of reads to the reference sequence [10]\n";
    print "-ws2 INT     minumun word size alignment of reads to the LTR sequence [10]\n";
    print "-bp          path to alternative blastn installation\n";
    print "-kt          keep reads database files (in case of very large input reads file) \n";
    
    exit;
}


if ($uref)
{
    print "I will use user provided fasta file as reference.\n"
}
elsif ($coo)
    {
        print "Referece sequence will be extracted from coordinates ".$coo."\n";
        print "from raw fasta reference in folder: ".$hg19_ref."\n";
        print "LTR referece sequence in folder: ".$vir_ref."\n";
        print "Reads to align to reference: ".$reads."\n";
    }
else
{
  print "Please provide suitable input files. Check usage with breakalign_5_5_hg19_hiv.pl --help\n";  
}

#
#exit;

#candidate region genome coordinates, input must be like chr3:76076289-76076636
#my $coord_ref = $ARGV[0];
my $fasta_ref;

unless ($uref)
{
    my @ar_c = split(":",$coo);
    #die $coo;
    my $chr_c = $ar_c[0];
    my @st_en = split("-",$ar_c[1]);
    my $st = $st_en[0];
    my $en = $st_en[1];


    $fasta_ref = "hg19_".$chr_c."_".$st."_".$en.".fa";
    #die $fasta_ref;
    open(FASTAREF,">$fasta_ref") or die "I cannot open $fasta_ref";
    print FASTAREF ">hg19_".$chr_c."_".$st."_".$en."\n";
    my $seqio = Bio::SeqIO->new(-file => $hg19_ref.$chr_c.".fa", -format => "fasta");
    #my $seqio = Bio::SeqIO->new(-file => "/databank/raw/hg19_full/".$chr_c.".fa", -format => "fasta");
    while(my $seq = $seqio->next_seq)
    {
        #my $refer = $seq->subseq($st,$en);
        #uc $refer;
        print FASTAREF uc($seq->subseq($st,$en)), "\n";
        print uc($seq->subseq($st,$en)), "\n";
    }

    close (FASTAREF);
}
else
{
    $fasta_ref = basename($uref);
    #copy($uref,$input_dir);
}


#die;

#fasta file with reads sequences, possibly not all of them, only the ones mapping the candidate region or knowing to contain the LTR
my $reads_f = basename($reads);
#my $reads_f = $ARGV[1];
my $word_size = $word_s;
#my $word_size = $ARGV[2];
my $word_size_b = $word_s2;
#my $word_size_b = $ARGV[3];
#unless ($ARGV[2]) {$word_size = 10}; #16/08/13 changed word size interactive even for the first blast 


#unless ($ARGV[3]) {$word_size_b = 20};

#label for unique output
my $label;
unless ($uref)
{
    $label = $fasta_ref;
    $label =~  s/\.[^.]+$//;
    #$label =~ s/.fa$//;
}
else {$label ="User_reference"}
#die $label;

#Generating a new input file with unique and numeric ids

my $new_reads_f= "Unique_".$reads_f;

open(FH, '>', $new_reads_f) or die "Could not open file '$new_reads_f' $!";

open(FIN, $reads) or die "Could not open file '$reads' $!";

my $n_line = 0;
while (<FIN>)
{
    my $line = $_;
    chomp($line);
    unless($line) {next;}
    
    $line =~ s/\r\n$//; #just in case chomp would not work
    $line =~ s/\r$//;#just in case...to deal with any os
    if ($line =~ /^>/)
    {
        $n_line++;
        $line = ">".$n_line;
        
    }
    print FH $line."\n";
}


$reads_f = $new_reads_f;


#was blast path on deva
#my $blast_path = "/package/blast/2.2.23/bin/";

#blast path to a more recent version
#my $blast_path = qx(which blastn);
#my $blast_path = "/package/blast/ncbi-blast-2.2.29+/bin/";

#die $blast_path;

#blast default on cruncher, better not to use it now
#my $blast_path = "/usr/bin/";

#input and output folder
my $input_dir = getcwd;
my $output_dir = getcwd;

if ($input_dir =~ m/ /) {$input_dir =~ s/ /\\ /};
if ($output_dir =~ m/ /) {$output_dir =~ s/ /\\ /};

#copy the reference sequence in the input dir in case it is not already in the working directory
copy($uref,$input_dir);
#to handle escape character in -e or open in perl
my $output_dir_b =  $output_dir;
if ($output_dir_b  =~ m/\\ /) {$output_dir_b =~ s/\\ / /;}


# my temp output
my $temp_out = $output_dir."/temp_br_al";
#my report output
my $rep_f = $output_dir."/Report_".time().".txt";
# I will label temporary the main output with a time stamp

#my $out_file = $output_dir."/align_for_breakpoint_".time()."_".$label;
my $out_file = $output_dir."/align_for_breakpoint"."_".$label.".txt";
#blast output file
my $blast_f = $temp_out."/blast_out_ref_";
#breakalign fasta input
#my $br_input =$output_dir."/hits_seq_".$label.".fa";
my $br_input =$temp_out."/hits_seq_".$label.".fa";
#my $br_input =$temp_out."/hits_seq_".time()."_";
my $br_input_2 =$temp_out."/hits_seq_K03455_".$label.".fa";
#my $br_input_2 =$temp_out."/hits_seq_K113_".time()."_";

#open output files

if ($rep_f =~ m/\\ /) {$rep_f =~ s/\\ / /};
open FILER, ">$rep_f" or die "unable to open $rep_f $!";

if ($out_file =~ m/\\ /) {$out_file =~ s/\\ / /};
open(OUT_F,">$out_file") or die "I cannot open $out_file\r\n" ;

#create a temp folder if there is not one
if ($temp_out =~ m/\\ /) {$temp_out =~ s/\\ / /}; # because -e does not want escape character
unless (-e $temp_out)
{
    my $cmd_1 = "/bin/mkdir $temp_out";
    #die print $current_path;
    print "executing: $cmd_1\n";
    print FILER "executing: $cmd_1\n";
    system($cmd_1) == 0 or die "system cmd [$cmd_1] failed ($?): $!";          # fork, execute $cmd and return exit status
}

#Running makeblastdb to make the blast database of the candidate region in the correct format 
#EX: makeblastdb -in ./pre-int_K_15.fa -dbtype nucl -out ./pre_int_K_15 -title "pre_int_K_15"

my $input_dir_b = $input_dir;

if ($input_dir =~ m/\\ /) {$input_dir_b =~ s/\\ / /;};
unless (-e $input_dir_b."/blastdb_reads.nin" or -e $input_dir_b."/blastdb_reads.00.nin" or -e $input_dir_b."/All_reads_db.00.nhr")
{

    if (-e $input_dir_b."/".$reads_f)
    {
        
        ###here i am going to make the database using the reads, not the genomic region
        #my $cmd_2 = $blast_path."/makeblastdb -in ./".$reads_f." -dbtype nucl -parse_seqids -out ./All_reads_db -title 'All_reads_db'";
        my $cmd_2 = $blast_path."/makeblastdb -in ".$input_dir."/".$reads_f." -dbtype nucl -parse_seqids -out ".$output_dir."/All_reads_db -title 'All_reads_db'";
        
        print "executing: $cmd_2\n";
        print FILER "executing: $cmd_2\n";
        system($cmd_2) == 0 or die "system cmd [$cmd_2] failed ($?): $!";
            
    }
    else
    {
        print "The reference file $reads_f does not exist or it is not in the input directory: ".$input_dir;
    }
}
else {print "Blast database blastdb_reads exists in the current directory, it will not be redone\n"}

#Running blast using the reads sequences in input versus the candidate region
#Ex: blastn -query ../breakalign_tests/kanapin_cand.fas -db "pre_int_K_15" -word_size 20 -outfmt 4 -out ../breakalign_tests/blast_out_pre_int_K_15_all

# i was here make the check for doing blastn with All_reads_db!!!!!!~~~~~~~~~~~

my $blast_f_label =  $blast_f.$label;
if ( $blast_f_label =~ m/\\ /) {$blast_f_label =~ s/\\ / /};

unless (-e $blast_f_label)
{
    my $cmd_3 = $blast_path."/blastn -query ".$input_dir."/".$fasta_ref." -num_threads 8 -db 'All_reads_db' -word_size $word_size -outfmt '7 sseqid qstart qend sstart send qseq sseq evalue bitscore score length qframe sframe' -out ".$blast_f_label; 
    print "executing: $cmd_3\n";
    print FILER "executing: $cmd_3\n";
    system($cmd_3) == 0 or die "system cmd [$cmd_3] failed ($?): $!";

}
else
{
    print "The genomic region in file $label has been analyzed with blastn before."."\n"."The output file is in ".$blast_f_label."\n"; 
    #die "The genomic region file $fasta_ref does not exist or it is not in the input directory: ".$input_dir;
}
######################################################################################################
##I am converting the reference into an array and then into a hash, I would be able to handle indels##
######################################################################################################

open(FA, $fasta_ref) or die "I cannot open $fasta_ref";
my @array_ref;
my %hash_ref;
my $ref;
while(<FA>)
{
    my $line = $_;
    chomp($line);
    $line =~ s/\r\n$//; #just in case chomp would not work
    $line =~ s/\r$//;#just in case...
    if ($line =~ /^>/) {next};
    $ref .= $line; #it should handle new lines characters in the reference sequence 
}

close(FA);

@array_ref = split('',$ref);

my $i =1 ;
foreach (@array_ref)
{
    #maybe just simply $hash_ref{$i}= $_;
    $hash_ref{$i}{$_} = ''; #so each position can be either a base or a gap '-'
    $i++;
    
};
######################################################################################################

#####################################
#Storing reads with hits into a hash# --> hits in the region or the fasta file provided
#####################################
 
# I won't scann for the significant read resequeneces in the file in input, it could be too big, I'll use blastdbcmd instead to retrieve the sequence with hits
#now checking what read ids have singinificant hits and keep the sequence
my %hits;
my %query_al;
my %read_al;
#opening fasta file that will be the input of breakalign
#open(HITS_F,">>$br_input") or die "I cannot open $blast_f\n";

#$blast_f_label = $blast_f.$label;

if ($blast_f_label =~ m/\\ /) {$blast_f_label =~ s/\\ / /};
my %seen;
open(BLASTF,$blast_f_label) or die "I cannot open $blast_f_label\n";
while (<BLASTF>)
{
  #die "I am here";
    if ($_ =~ /^#/) {next;}
    my $line = $_;
    chomp($line);
    my @ar = split('\t',$line);
    #this is to get only the best hit of each alignment
    if (exists($seen{$ar[0]})){next;}
    
    if ($ar[3] < $ar[4])
    {
        $hits{$ar[0]}{'sense'} = 'F'; #the alignment is forward
        $hits{$ar[0]}{'st_al'} = $ar[1];
        $hits{$ar[0]}{'st_al_read'} = $ar[3];
        #my $rsir = ($ar[1] - $ar[3])  + 1; #read start in reference
        #$hits{$ar[0]}{'rsir'} = $rsir;
        
    }
    else
    {
        $hits{$ar[0]}{'sense'} = 'R'; #the alignment is with the reverse complement of the read
        $hits{$ar[0]}{'st_al'} = $ar[1];
        $hits{$ar[0]}{'end_al_read'} = $ar[3];
       
    }
    my $align_q = $ar[5];
    my @ar2 = split('',$align_q);
    my $p = $ar[1]; #alignent on the query sequence
    foreach (@ar2)
    {
           my $base = $_;
           $query_al{$ar[0]}{$p}{$base} = '';
           $p++;
    }
    #if ($ar[0] eq '35'){die Dumper(\%query_al)}; "hyphone, '-' is stored
    
    my $align_r = $ar[6];
    my @ar3 = split('',$align_r);
    my $p2 = $ar[1]; #I put the read alignment in the query position coordinates
    foreach (@ar3)

    {
       my $base = $_;
       $read_al{$ar[0]}{$p2}{$base} = '';
       $p2++;
    }
    
    #now for each read id we get the sequence using blastdbcmd
    #unless ( -e "/hts/data2/emarchi/Retroseq_out/OV_participants/p_id_206ce0ed/temp_br_al/hits_seq_1354215609_hg18_chr17_4899617_4900052_fromR.fa")
    #unless ( -s $br_input)
    #{
        
        my $cmd_4 = $blast_path."/blastdbcmd -entry 'lcl|".$ar[0]."' -db 'All_reads_db' >> ".$br_input;
        #die $cmd_4;
        print "executing: $cmd_4\n";
        print FILER "executing: $cmd_4\n";
        system($cmd_4) == 0 or die "system cmd [$cmd_4] failed ($?): $!";
    #}
    $seen{$ar[0]} = '';
    
   
}
close(BLASTF);

#die Dumper(\%hits); #ok!

#1. Getting the sequenqueces with significant hits
#2. Blasting the sequences with hits versus K03455 (making a db of them)
#3. Producing a second filtered input for breakalign with only reads having hits with K03455, likely to contain the LTR
#IMPORTANT: WE HAVE TO MAKE SURE THAT EACH SEQUENCES HAS A UNIQUE ID
#IMPORTANT DIFFERENCE WITH PREVIOUS VERSION, HERE I OPEN THE FASTA FILE WITH SEQUENCES RETRIVED  BY BLASTDBCMD!!!
#open(RE,$reads_f) or die "I cannot open $reads_f";
my $br_input_b =  $br_input;
if ($br_input_b  =~ m/\\ /) {$br_input_b =~ s/\\ / /;}
open(RE,"$br_input_b") or die "I cannot open "."$br_input_b";

#FIRST REMOVING THE NEW LINES, THE  SCRIPT BELOW EXPECT FASTA FORMAT WITH NO NEW LINES CHARACTERS

my $conc_line;
my $firstline = 1;
my $outvar;
my $IDD;
while(<RE>)
{
    my $line = $_;
    chomp($line);
    $_ =~ s/\n$//; # in case chomp would not work...
    $line =~ s/\r$//; #remember to remove \r too because we are in bloody mac!!!
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
    $conc_line .= $line; #it should handle new lines characters in the reference sequence 
    
}
$outvar.= $IDD."\n".$conc_line."\n"; #this is to print the last line
close(RE);
open(REWR, ">".$br_input_b) or die "I cannot write on NEW_".$br_input_b."\n";
print REWR $outvar;
close(REWR);

#############################################################
####finding the hits for k113 pr HIV in the mapping reads####
#############################################################

if (-e $br_input_b)
    {
        #die $br_input; 
        ###here i am going to make the database using the reads, not the genomic region
        my $cmd_5 = $blast_path."/makeblastdb -in ".$br_input_b." -dbtype nucl -parse_seqids -out ./hits_reads_db_".$label." -title 'hits_reads_db_".$label."'";
	#bloody spaces in the paht causing problems again, now makeblastdb doesn't like it
	# my $cmd_5 = $blast_path."/makeblastdb -in ".$br_input." -dbtype nucl -parse_seqids -out ".$output_dir."/hits_reads_db_".$label." -title 'hits_reads_db_".$label."'";
        print "executing: $cmd_5\n";
        print FILER "executing: $cmd_5\n";
        system($cmd_5) == 0 or die "system cmd [$cmd_5] failed ($?): $!";
            
    }
    else
    {
        print "The fasta file $br_input with read sequences having hits to the candidate regions does not exist";
    }

#my $output_dir_b =  $output_dir;
#if ($output_dir_b  =~ m/\\ /) {$output_dir_b =~ s/\\ / /;};)

if (-e $output_dir_b."/hits_reads_db_".$label.".nhr")
{
    #my $cmd_6 = $blast_path."/blastn -query /home/manu/HIV_ref/K03455_LEFT_LTR.fasta -num_threads 8 -db 'hits_reads_db_".$label."' -word_size ".$word_size_b." -outfmt '7 sseqid qstart qend sstart send qseq sseq evalue bitscore score length qframe sframe' -out ".$temp_out."/K03455_blast_out_".$label;
    my $cmd_6 = $blast_path."/blastn -query ".$vir_ref." -num_threads 8 -db 'hits_reads_db_".$label."' -word_size ".$word_size_b." -outfmt '7 sseqid qstart qend sstart send qseq sseq evalue bitscore score length qframe sframe' -out ".$temp_out."/K03455_blast_out_".$label;
    print "executing: $cmd_6\n";
    print FILER "executing: $cmd_6\n";
    system($cmd_6) == 0 or die "system cmd [$cmd_6] failed ($?): $!";

}
else
{
    die "The second database of reads mapping the genomic region ".$output_dir_b."/hits_reads_db_".$label." file does not exist or it is not in the correct directory: ".$output_dir_b;
}
#generating the second breakalign input with only reads having hits toward K03455

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
    print "executing: $cmd_7\n";
    print FILER "executing: $cmd_7\n";
    system($cmd_7) == 0 or die "system cmd [$cmd_7] failed ($?): $!";
    #sono qui
}
close(BLASTH);

#FIRST REMOVING THE NEW LINES (as I did before), THE  SCRIPT BELOW EXPECT FASTA FORMAT WITH NO NEW LINES CHARACTERS

my $br_input_2_b =  $br_input_2;
if ($br_input_2_b  =~ m/\\ /) {$br_input_2_b =~ s/\\ / /;}


if (-e $br_input_2_b)
{
    open(RE2,"$br_input_2_b") or die "I cannot open "."$br_input_2_b";
    #open(RE2,"/hts/data2/emarchi/Retroseq_out/OV_participants/p_id_206ce0ed/temp_br_al/hits_seq_1354215609_hg18_chr17_4899617_4900052_fromR.fa") or die "I cannot open "."/hts/data2/emarchi/Retroseq_out/OV_participants/p_id_206ce0ed/temp_br_al/hits_seq_1354215609_hg18_chr17_4899617_4900052_fromR.fa";
    my $conc_line_2;
    my $firstline_2 = 1;
    my $outvar_2;
    my $IDD_2;
    while(<RE2>)
    {
        my $line = $_;
        chomp($line);
        $_ =~ s/\n$//; # in case chomp would not work...
        $line =~ s/\r$//; #remember to remove \r too because we are in bloody mac!!!
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
        $conc_line_2 .= $line; #it should handle new lines characters in the reference sequence 
        
    }
    $outvar_2.= $IDD_2."\n".$conc_line_2."\n"; # this is to print the last line
    close(RE2);
    open(REWR2, ">"."$br_input_2_b") or die "I cannot write on NEW "."$br_input_2_b"."\n";
    print REWR2 $outvar_2;
    close(REWR2);
    #### Second input to breakaling done
}
##########################################################################
#getting read ids and sequences and make the reverse complement sequences# ---> in this version I keep only the reads having hits to the virus sequence
##########################################################################


open(RE,$br_input_2_b) or die "I cannot open $br_input_2_b\n"; # $br_input_2_b contain reads found to have hits to the region and the virus at the same time
my $r_id;
#my %reads;
my %h_reads_s; #hashes with the read sequences saved position by position, like the reference
my %h_reads_rc; # i don't really need that and I'm not using it...but i might use it to solve the problem og bugs for too long printing...
    
my $length_re; #important to work out the proper alignement
while(<RE>)
{
    my $line = $_;
    chomp($line);
    $_ =~ s/\n$//; # in case chomp would not work...
    $line =~ s/\r$//; #remember to remove \r too because we are in bloody mac!!!
    $line =~ s/\r\n$//; # not needed but one or the other should work
    if ($line =~ /^>/)
        {
            #print $line."\n";
            if ($line =~ /(\d+)\s/)
            {
            #$line =~ s/^>lcl\|\d+ //; #last modified 8/11/16, this line was important but not anymore?if I do check about blastn versions
            $r_id = $1;
	    #$r_id = $line;
            
            #$r_id =~ s/^>//;
            #$r_id =~ s/^>lcl\|//;
            #die $r_id; # it was sticking an lcl| in each header 08/11/16
            }
            
            #$r_id =~ s/\#0$//; #important to remove the '#' from the key!!! another difference with perl using MAC!?
        }
        else
        {
            my $seq = $line;
            #die Dumper(\%hits);
            if (exists ($hits{$r_id})) #storing only read sequences having significant hits
            {
                #die "I get here";
                $length_re = 0; #getting the length of the read sequence
                my @ar = split('',$seq);
                foreach (@ar)
                {
                    $length_re++;     
                }
                            
                my $rsir;
                if ($hits{$r_id}{'sense'} eq 'F')
                {
                    my $rsir = ($hits{$r_id}{'st_al'} - $hits{$r_id}{'st_al_read'})  + 1; #read start in reference
                    my $i=$rsir;
                    $length_re = 1; #just to work out the length of the read
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
                    #my $rsir = $hits{$r_id}{'st_al'} - ( $hits{$r_id}{'end_al_read'} - $length_re);
                    my $ii=$rsir;
                    foreach (@ar_rc)
                    {
                        $h_reads_s{$r_id}{$revc}{$ii}{$_}= $_;
                        
                        #$h_reads_rc{$r_id}{$seq}{$ii}{$_}= $_;
                        $ii++;
                        #$ii++;
                    }
                }
                
                
 
                
                
            }   
            
        }
    
    }

close(RE);

#die Dumper(\%h_reads_s); # ok
#########################################################
############Printing out the alignment ##################
#########################################################

#the genomic region first
###making the header with the reference sequence
foreach my $pos ( sort _numerically_asc keys %hash_ref)
{
    
    foreach my $b (keys %{$hash_ref{$pos}})
    {
        print OUT_F $b;
        
    }
    
    
}

print OUT_F "\n"; 

#die Dumper(\%h_reads_s);

my $out; # output variable to print out
#my $line_out;
my $indent;
foreach my $ids (sort keys %hits)
{
    #if (exists($seen{$ar[0]})){next;}
    unless (exists($h_reads_s{$ids})){next;}  # -> this is just to print out the reads found to have hits

    foreach my $seq (keys %{$h_reads_s{$ids}}) #just getting the indentation using the query alignment start, whish is the start of the alignment on the genomic region I'm testing
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
    foreach my $seq (keys %{$h_reads_s{$ids}}) #checking here which one of the four possibile sequences for each id to print (the mate that align and in which direction)
    {
        
        foreach my $pos (sort _numerically_asc keys %{$h_reads_s{$ids}{$seq}}) 
        {
            #unless (exists ($read_al{$ids}{$pos})) #changed for read 35 bug with insertion 22/05/20
            unless (exists ($query_al{$ids}{$pos}))
            {
                foreach my $b (keys %{$h_reads_s{$ids}{$seq}{$pos}}) #taking the only base possible there
                {
                    $line_out .= lc($b);
                }              
            }
            else
            {
                #foreach my $b (keys %{$read_al{$ids}{$pos}}) #taking the only base possible there
                foreach my $b (keys %{$query_al{$ids}{$pos}})#taking the only base possible there ##changed for read 35 bug with insertion 22/05/20
                {
                    #unless($b eq '-'){$out .= uc($b)};
                    if($b eq '-')
                    {
                      $line_out =~ s/^\s//;#just removing a white space from the indentation
                      $inser.=(keys %{$read_al{$ids}{$pos}})[0]; #in the same position I check what base should be
                      $b = (keys %{$read_al{$ids}{$pos}})[0]; #22/06/20 changed and print the base instead the "-" as in read alignment as Robert prefer, not sure is ideal
                      $line_out .= lc($b); #I would mark and make visible the insertion with a lower case letter
                      next;# this is just to the next base otherise a line below would make an upper case again
                    } 
                    if((keys %{$read_al{$ids}{$pos}})[0] eq '-') #checking if there is deletion and it would be showed in the read alignment with a "-"
                    {
                      $dele.=(keys %{$query_al{$ids}{$pos}})[0]; #in the same position I check what base should be
                      
                      $b = (keys %{$read_al{$ids}{$pos}})[0] #22/06/20 same as above change, but I would print the "-" instead of the base (from query alignent anyway) because in the read is missing anyway
                    }
                    $line_out .= uc($b); # the "-" character is not affected
                }        
            }
            
        }
   
    }
    #if ($inser ne ''){ $line_out .= " Insertion detected: ".$inser." "; $inser = '';}
    #if ($dele ne ''){$line_out .= " Deletion detected: ".$dele." ";$dele = ''}
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
#Delete temporary files#
########################


unless ($k_f)
{
  #remove temp files but in version 5.6 I removed temp directory with hits anyway, consider that...
  opendir(DIR, $temp_out);
  my @t_files = readdir DIR;
  
  foreach my $t_file (@t_files)
  {
   unlink($temp_out."/".$t_file);
   print "Deleteted temporary file: $t_file\n";
  }
  closedir(DIR);
  rmdir($temp_out)  or warn "couldn't remove $temp_out: $!";

  opendir(DIR, $input_dir);
  my @files = grep(/^All_reads_|^hits_reads_/,readdir(DIR));
  foreach my $file (@files)
  {
   unlink($input_dir."/".$file);
   print "Deleteted reads database file: $file\n";
  }
  closedir(DIR);
}

#############################################################################################


#making the html output reusing  br_text_to_html_mark.pl

br_txt_to_html($out_file);

#####################################
##########Subroutines################
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
        
        if ($line =~ /^[A-Z]/) #must be only the first line, the reference
        {
            $html_out.= $line."<br>";
            next;
        }
        
        if ($line =~ /(^\s+)/)
        {
           
            $indent = $1;
            $line =~ s/^\s+//;#code
            
            if ($line =~ /(\s+\d+)/) #the id are the only numeric and not the last
            {
                $r_id = $1;
                $line =~ s/\s+\d+//;
                #if ($line =~ m/(\*\s[INSERTION:|DELETION:]\s[A-Z]+)/)
                if ($line =~ m/(\*\s.+$)/)
                {
                  #die $r_id;
                  $r_id .= $1;
                  $line =~ s/\*.+$//;
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




