
use warnings;
use strict;


use Bio::SeqIO::fasta;
use Bio::SearchIO;

my $inputData=$ARGV[0];
my $outputdir=$ARGV[1];
my $cutoffLength=$ARGV[2];
my $prefix=$ARGV[3];
my $numArgs = $#ARGV + 1;

print "$numArgs,  Output directory is $outputdir, input dataset is $inputData, cutoff length is $cutoffLength \n";

open(IN,  "$inputData") or die "cannot open input file\n";

my $resFile2 = ${outputdir}."/".$prefix."_grt".$cutoffLength.".fa";
open (OUT3, " > $resFile2") or die "cannot open output file\n";
print "FASTA Output in $resFile2\n";

my $seqio = Bio::SeqIO->new(-file=> $inputData,-format => 'Fasta' );

my $sum_contig_length=0;
while ( my $seq_obj = $seqio->next_seq()){
    my $read_id=$seq_obj->id;
    my $read_seq=$seq_obj->seq;
     
    
    my $realLength = (split /\_/, $read_id)[-5];
    my $contigID = (split /\_/, $read_id)[-7];

    

    if ($realLength>=$cutoffLength){
	print OUT3 ">",$read_id, "\n", $read_seq, "\n";
    }
}





