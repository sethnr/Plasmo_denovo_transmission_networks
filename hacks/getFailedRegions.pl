#!/usr/bin/perl

while(<>) {
    chomp;
    unless ($_ =~ m/^\#\d+\:DONE/) {
	@F = split;
	print $F[-1]."\n";
    }
}


