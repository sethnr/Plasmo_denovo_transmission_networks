#!/usr/bin/perl

$targetFile = shift;
$strFile = shift;

open(STRS,$strfile);
open(TARGETS,$targetFile);

while(<TARGETS>) {
	      my $name, $chr, $st, $en = split(/\s/,$_);
	      push @{$targets{$chr}}, \($name,$st,$en);
}

print %targets;

while(<STRS>) {
  my @F = split(/\s/,$_);
  foreach $target @{$targets{$F[0]}} {
    ($name, $st, $en) = @{$target};
  }
}
