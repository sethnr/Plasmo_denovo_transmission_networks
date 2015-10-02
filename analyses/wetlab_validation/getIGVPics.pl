#!/usr/bin/perl

my ($REF, $BED, @SAM) = @ARGV;

my $IGV="~/software/IGV_2.3.59/igv.sh";
my $IGVOUT=$BED;
$IGVOUT =~ s/\.bed$/.igv/gi;

open(BED,$BED);
open(IGV,">$IGVOUT");

print IGV "genome $REF\n";
foreach my $SAM (@SAM) {
    print IGV "load $SAM \n";
}

while(<BED>) {
    print $_;
    my ($chr,$st,$en,$name) = split("\t",$_);
    next if $chr eq '@SQ';
    $name =~ s/\s//gi;
    $name =~ s/\W/_/gi;

    print IGV <<"IGVEND";
goto $chr:$st-$en
sort base
collapse
snapshot $name.png

IGVEND

}
