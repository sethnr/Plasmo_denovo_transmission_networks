#!/usr/bin/perl

my @lastindel = ();
my $lastpos = -1;
while(<>) {
    next if ($_ =~ m/^\#/);
    my @F = split;
#    print $F[1].", 1\n";
#    print "@F[0..4] , 1-4\n";
    my ($chrom, $pos, $id, $ref, $alt) = @F[0..4];
    $pos = $pos+0;
    if (length($ref) > 1 || length($alt) > 1) {
#	print "indel",$pos,":",$lastpos."\n";
	if (($pos == $lastpos+1) || ($pos == $lastpos)) {
	    print "neighbour\n";
	    print join("\t",@lastindel)."\n";
	    print join("\t",@F)."\n";
	}
	@lastindel = @F;
	$lastpos = $pos;
    }
}
