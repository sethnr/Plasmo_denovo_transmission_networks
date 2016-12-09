#!/usr/bin/perl

import re;



while(<>) {
    my @F = split(/\s+/,$_);
    $cigar = $F[5];
    
    while ($cigar =~ s/^(\d+)([MSID])//){
#	print $1."\t".$2."\t".$cigar."\n";
	$c{$2} += $1;
	if ($2 eq "S" && $1 > 2) {} else {$c{"ALL"} += $1;}
    }
}
print STDOUT $ARGV."\t".
    $c{"M"}."\t".$c{"I"}."\t".$c{"D"}."\t".$c{"S"}.
    "\t".$c{"ALL"}."\n";
