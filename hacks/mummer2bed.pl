#!/usr/bin/perl/

my $coords = shift;
open(COORDS,$coords);

while(<COORDS>){
    my @F = split;
    if ($#F==1) {
	
	($F1, $F2) = @F;
	print STDERR $F1."\t".$F2."\n";

	$F1 =~ s/\.fa$//i;
	$F1 =~ s/\.fasta$//i;
	$F1 =~ s/^.*\///i;
	$F2 =~ s/\.fa$//i;
	$F2 =~ s/\.fasta$//i;
	$F2 =~ s/^.*\///i;
	$fn1 = "nucmer.".$F1."in".$F2.".bed";
	$fn2 = "nucmer.".$F2."in".$F1.".bed";
	print STDERR $fn1."\t".$fn2."\n";

	open(F1INF2,">$fn1");
	open(F2INF1,">$fn2");
	print F1INF2 join("\t","chrom","chromStart","chromEnd")."\n";
	print F2INF1 join("\t","chrom","chromStart","chromEnd")."\n";
    }
    elsif ($F[0] =~ m/^\=+/gi) {next;}
    elsif ($#F==12) {
	for (my $i=0; $i <= $#F; $i++) {
	    $F[$i] =~ s/\s//gi;
	}
	my ($S1,$E1, undef, $S2, $E2, undef, $L1 ,$L2, undef , $PCID, undef, $C1, $C2) = @F;
	print F1INF2 join("\t",$C1,$S1,$E1)."\n";
	print F2INF1 join("\t",$C2,$S2,$E2)."\n";
    }
}
