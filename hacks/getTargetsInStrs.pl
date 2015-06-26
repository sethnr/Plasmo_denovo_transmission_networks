#!/usr/bin/perl

$targetFile = shift;
$strFile = shift;

open(STRS,$strFile);
open(TARGETS,$targetFile);

while(<TARGETS>) {
#  print $_;
#  my ($name, $chr, $st, $en) = split(/\s+/,$_);
  my ($chr, $st) = split(/\s+/,$_);
#  my @target = ($name,$st,$en);
  my @target = ($name,$st,$st);
#  print "@target";
#  print "$name, $chr, $st, $en !\n";
  push @{$targets{lc($chr)}}, \@target;
}

# foreach $chr (keys %targets) {
#   foreach $target (@{$targets{$chr}}) {
#     print "\#".join("\t", $chr, @{$target})."\n";
#   }
# }

while(<STRS>) {
  my @F = split(/\s+/,$_);
#  print join("\t",$F[0], $F[1], $F[3],$_);
  foreach $target (@{$targets{lc($F[0])}}) {
    ($name, $st, $en) = @{$target};
    if ($F[1] < $en && $F[3] > $st) {
      print "\n\#\#$name:$st-$en\n" unless $targeted{$name};
      $targeted{$name} = 1;

      print $_;
    }
  }
}
