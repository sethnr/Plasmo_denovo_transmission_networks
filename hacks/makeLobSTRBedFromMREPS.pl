#!/usr/bin/perl 

my $copy_number = 1;

while(<>) {
    my @F = split;
    my ($chrom, $start, undef, $end, undef, $size, $period, $exponent, $error, $repunit) = @F[0..9];
    if(length($repunit) != $period) { $repunit = $F[10];}
    $period =~ s/\D//gi;
    $exponent =~ s/[\[\]]//gi;
    $exponent = $exponent +0;
    $period = $period +0;
    $start = $start +0;
    $end= $end+0;
    $length = $end-$start;
    print join("\t",
		   $chrom,
	       $start,
	       $end,
	       $period,
	       $copy_number,
	       ".",".",".",
#	       $period*$exponent,
	       $length*2,
	       ".",".",".",".",".",
	       $repunit)."\n";
}
