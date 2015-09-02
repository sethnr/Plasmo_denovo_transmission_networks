#!/usr/bin/perl

use POSIX strftime;

my $date = strftime '%y%m%d', localtime;

my $outfile= $ARGV[0];

$outfile =~ s/sub_/done_/;
$outfile =~ s/out.log/$date.txt/;
print $outfile."\n";

print $date."\n";

open(OUT,">$outfile");
while (<>) {
    my @F = split(/\s+/,$_);
    if ($F[0] eq '#DONE:'){
	print OUT join(" ",@F[1..$#F])."\n";
    }
}

exit 0
