#!/usr/bin/perl

my %combs;
my %inds;

while(<>) {
    next if $_ =~ m/^#/;
    @F = split;
    $filters = $F[6];
    
    $combs{$filters}++;
    my @filters = split(";",$filters);
    foreach $f (@filters) {
	$inds{$f}++;
    }
}

print "INDIVIDUAL FILTERS\n";
foreach $f (sort {$inds{$b} <=> $inds{$a}} keys %inds) {
    print "\t".$f."\t".$inds{$f}."\n";
}
print "\n\n";
print "COMBINED FILTERS\n";
#foreach $f (sort {$combs{$b} <=> $combs{$a}} keys %combs) {
foreach $f (sort {$a cmp $b} keys %combs) {
    print "\t".$f."\t".$combs{$f}."\n";
}
