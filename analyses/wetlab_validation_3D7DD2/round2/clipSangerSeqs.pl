#!/usr/bin/perl

#my %clip, %cut;


my $primerfile = shift;

open(PRIMERS,$primerfile);

while (<PRIMERS>) {
  print $_;
  ($id, $chr, $var, $len, $F, $R) = split(/\s+/,$_);
  $clip{$id."F"}=$F;
  $cut{$id."F"}=$len;
  print join("\t",$id,$F,$R)."\n";
  $R =~ tr /atcgATCG/tagcTAGC/; $R = reverse($R);
  $clip{$id."R"}=$R;
  $cut{$id."R"}=$len;
  print join("\t",$id,$F,$R)."\n";
    
}


foreach $PRIM (keys(%clip)) { 
#    print $PRIM;
    @files = `ls *${PRIM}.seq`;
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
