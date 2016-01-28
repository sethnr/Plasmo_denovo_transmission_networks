#!/usr/bin/perl

#my %clip, %cut;

$clip{"1F"}="GgaaaaaaAACAAGAGGGGGAGT";
$clip{"1R"}="CCCCCTCGACGATGTAGTTG";
$clip{"2F"}="GGACATGTCTCCCAAGGTCG";
$clip{"2R"}="AACGCTTTGACCTTAAAAGATGT";
$clip{"3F"}="AAGCACGTTCCCTTAAGTGT";
$clip{"3R"}="TGCGACCTGAAATTACACCTCA";
$clip{"4F"}="TTGGGGTACAGGGAAATCTATATTT";
$clip{"4R"}="TGATGCAAAAAAGGGATTGAATTTG";
$clip{"5F"}="TGTGGGGGTCCTTATGATCT";
$clip{"5R"}="ACCCCCTTTGGATGGTTCAG";
$clip{"6F"}="CATTCCTCCCACCAATGGCT";
$clip{"6R"}="TGTTGTAAGGCAGTGTGTTGC";
$clip{"7F"}="CGTATGTTCGTCTTCTTCATTCCA";
$clip{"7R"}="TGCAAAGAATTCGAAAGTTCAGA";
$clip{"8F"}="ATTCAATTGCTGCGTGACCC";
$clip{"8R"}="AGAACATCGCGATTTATGTACTTCT";
$clip{"9F"}="TTCACGCCAAGTTTCATCAAC";
$clip{"9R"}="TGCATGTATATGGTTCTTTCCCG";
$clip{"10F"}="TCCCGAAGAAGAAGTTGTCACA";
$clip{"10R"}="GTGTTTACACAACTGTGGATGAGA";
$cut{"1F"}=894;
$cut{"2F"}=286;
$cut{"3F"}=789;
$cut{"4F"}=351;
$cut{"5F"}=1006;
$cut{"6F"}=223;
$cut{"7F"}=173;
$cut{"8F"}=315;
$cut{"9F"}=221;
$cut{"10F"}=197;
$cut{"1R"}=894;
$cut{"2R"}=286;
$cut{"3R"}=789;
$cut{"4R"}=351;
$cut{"5R"}=1006;
$cut{"6R"}=223;
$cut{"7R"}=173;
$cut{"8R"}=315;
$cut{"9R"}=221;
$cut{"10R"}=197;

foreach $PRIM (keys(%clip)) { 
#    print $PRIM;
    @files = `ls *${PRIM}.ab1.seq`;
#    print @files;
    foreach $FILE (@files) {
	chomp($FILE);
#        print STDERR $FILE;
#	print STDERR "fasta_formatter -i $FILE | fastx_clipper -d1 -a $clip{$PRIM}\n";
        $unclipped=`fasta_formatter -i $FILE`;
        $clipped=`fasta_formatter -i $FILE | fastx_clipper -d1 -a $clip{$PRIM}`;
	if (length($unclipped) == length($clipped)) {
#		print $unclipped;
#		print "TRUNCATED\t".$cut{$PRIM}."\n";
		@lines = split("\n",$unclipped);
#		print $lines[1];
		$lines[1] = substr($lines[1],0,$cut{$PRIM});
		print $lines[0]." TRUNCATED\n";
		print $lines[1]."\n";
		
#		print join("\n",@lines)."\n";
	    }
	else {
		#print $unclipped;
#		print "CLIPPED\n";
		@lines = split("\n",$clipped);
		print $lines[0]." CLIPPED\n";
		print $lines[1]."\n";
#		print $clipped;
	    }
	print "\n";
    }
}
