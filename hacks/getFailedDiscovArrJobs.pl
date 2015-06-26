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

open(SUBMITTED,"${FOLDER}/sub_commands.txt");
print "${FOLDER}/sub_commands.txt";

while(<SUBMITTED>) {
    my @F = split;
    $job = $F[-3]."_".$F[-2]."_".$F[-1];
    if (m/^\#\d+:DONE/gi) {
    print STDERR join("\t",$job,
		      1,$started{$job},$finished{$job} || '-')."\n";
#	print $_ unless $finished{$job}
    }
    else {
    print STDERR join("\t",$job,
		      1,$started{$job},$finished{$job})."\n";
#	print $_ unless $finished{$job}
    }
}
