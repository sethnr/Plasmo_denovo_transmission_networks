#!/usr/bin/perl

my $FOLDER = shift;

open(BAMFILES,"ls -1 ${FOLDER}/*/*bam |");
while(<BAMFILES>) {
    chomp;
    my ($group,$sample,$job) = split('/',$_);
    $job =~ s/\..*$//gi;
#    print STDERR $_."\n".$job."!\n";
    $started{$job} = 1;
}
close(BAMFILES);

open(VCFFILES,"ls -1 ${FOLDER}/*/*vcf.gz |");
while(<VCFFILES>) {
    chomp;
    my ($group,$sample,$job) = split('/',$_);
    $job =~ s/\..*$//gi;
#    print STDERR $_."\n".$job."!\n";
    $finished{$job} = 1;
}
close(VCFFILES);

#exit(0);

open(SUBMITTED,'submitted.log');
while(<SUBMITTED>) {
    next if m/^\#/gi;
    if($_ =~ m/-J\s+(\S+)\s+/gi) {
	$job = $1;
	print STDERR join("\t",$job,
		   1,$started{$job},$finished{$job})."\n";
	print $_ unless $finished{$job}
    }
    else {
	print STDERR "job name not parsed for \n".$_."\n";}
    
}
