#! /usr/bin/perl -w
#
# File: reichenowi.allels.V2.pl
# Time-stamp: <16-Oct-2013 13:32:55 tdo>
# $Id: $
#
# Copyright (C) 2013 by Pathogene Group, Sanger Center
#
# Author: Thomas Otto
#
# Description:
#


my $coords=shift;
my $snpC  =shift;
my $snp   =shift;
my $ref   =shift;
my $result=shift;

if (!defined($result)) {
die "usage: <filtered coords> <snp -C> <snp> <reference multifasta> <result>\n";

}

getMutations($snpC,$coords,$result);

my $ref_coords=getCoords($coords);
my $ref_snpC  =getSNPc($snpC);
my $ref_snp   =getSNP($snp);

my ($ref_seq,$ref_fai) = getFasta($ref);


foreach my $chr (keys %$ref_fai) {
  foreach my $pos (1..$$ref_fai{$chr}) {


	if (!defined($$ref_coords{$chr}{$pos})) {
	 # print "0\t0\t0\t-1";
	}
	else {
	  my $score=1;
	  my $base=substr($$ref_seq{$chr},($pos-1),1);
	  print "$chr\t$pos\t";

	  print "1\t$base\t";

	  ## confident calling
	  if (!defined($$ref_snpC{$chr}{$pos})) {
		print "$base\t";
	  }
	  else {
		print $$ref_snpC{$chr}{$pos}."\t";
		$score++;	  }


	  ### SNP call not so confident
	  if (!defined($$ref_snp{$chr}{$pos})) {
		print "$base";
	  }
	  else {
		print $$ref_snp{$chr}{$pos};
		$score--;
	  }
	  print "\t$score\n";

	}

#	print "\n";

  }
}

sub getCoords{
  my $file = shift;

  open (F, $file) or die "Problem to open SNP File: $file \n";

  my %h;

  ### get an field, where there is coverage
  foreach (<F>) {
	my ($refPos1,$refPos2,$queryPos1,$queryPos2,$overlap1,$overlap2,$identity,$refLength,$queryLength,$dum1,$dum2,$reference,$query) = split(/\s+/);

	foreach my $pos ($refPos1..$refPos2) {
	  $h{$reference}{$pos}=1
	}

  }

  close(F);
  return \%h;

}


sub getSNPc{
  my $file = shift;

  open (F, $file) or die "Problem to open SNP File: $file \n";

  my %h;
  foreach (<F>) {
	my ($refPos,$refWhat,$queryWhat,$queryPos,$dum1,$dum2,$refStrand,$queryStrand,$reference,$query) = split(/\s+/);

	###three cases
	if ($refWhat eq ".") {
	  ###insertion reference
	  if (defined($h{$reference}{$refPos})) {
		$h{$reference}{$refPos}.= "$queryWhat";

	  } else {
	  $h{$reference}{$refPos}.= "Ins: $queryWhat"
	}

	}
	elsif ($queryWhat eq ".") {
	  $h{$reference}{$refPos}="Del"
	}
	else {
	  $h{$reference}{$refPos}=$queryWhat
	}

  }

  close(F);
  return \%h;
}
sub getSNP{
  my $file = shift;

  open (F, $file) or die "Problem to open SNP File: $file \n";

  my %h;
  foreach (<F>) {
	my ($refPos,$refWhat,$queryWhat,$queryPos,$dum1,$dum2,$refStrand,$queryStrand,$dum3
,$dum4,$reference,$query) = split(/\s+/);

	###three cases
	if ($refWhat eq ".") {
	  ###insertion reference
	  if (defined($h{$reference}{$refPos})) {
		$h{$reference}{$refPos}.= "$queryWhat";

	  }else {
	  $h{$reference}{$refPos}.= "Ins: $queryWhat"
	}

	}
	elsif ($queryWhat eq ".") {
	  $h{$reference}{$refPos}="Del"
	}
	else {
	  $h{$reference}{$refPos}=$queryWhat
	}

  }

  close(F);
  return \%h;
}



sub getFasta{

  my $file = shift;

  open F, $file or die "Problmen\n";

  my %h;
  my %fai;

my $chr='';
my $length;

while ((<F>)) {
  chomp;

  if ((/>(\S+)/)) {
	$chr=$1;
  }
  else {
	chomp;
	$h{$chr}.=$_;
	$fai{$chr}+=length($_);
  }

}
return (\%h,\%fai);

}



sub getMutations{
  my $fileNameSNP     = shift;
  my $fileNameCoords  = shift;
  my $resultName      = shift;

  open (F, $fileNameSNP) or die "Problem to open SNP File: $fileNameSNP \n";

  my @File=<F>;
  close(F);

  my %BB;
  my %LB;

  my (%sizeRef,%sizeQuery,%coveredRef,%coveredQuery,%noCovRef,%noCovQry);
  my %h_sizeRef;
  my %h_sizeQuery;

  foreach (@File) {
	my ($refPos,$refWhat,$queryWhat,$queryPos,$dum1,$dum2,$refStrand,$queryStrand,$reference,$query) = split(/\s+/);

	if ($refWhat eq ".") {
	  $BB{$reference} .="$reference\tBBA\tIns\t$refPos\t$refPos\t0\t+\t.\tnote=\"Insert
ion+in+query:+$queryWhat++++Strand+of+query+is+$queryStrand\"\n";
	  $LB{$query} .="$query\tBBA\tDel\t$queryPos\t$queryPos\t0\t+\t.\tnote=\"Deletion+i
n+reference++++Strand+of+query+is+$queryStrand\"\n";
	}
	elsif($queryWhat eq "."){
	  $BB{$reference} .="$reference\tBBA\tDel\t$refPos\t$refPos\t0\t+\t.\tnote=\"Deleti
on+in+query++++Strand+of+query+is+$queryStrand\"\n";
	  $LB{$query} .="$query\tBBA\tIns\t$queryPos\t$queryPos\t0\t+\t.\tnote=\"Insertion+
in+reference:+$refWhat++++Strand+of+query+is+$queryStrand\"\n";
	}

	else {
	  $BB{$reference} .="$reference\tBBA\tSNP\t$refPos\t$refPos\t0\t+\t.\tnote=\"SNP+in
+query:+$queryWhat++++Strand+of+query+is+$queryStrand\"\n";
	  $LB{$query} .="$query\tBBA\tSNP\t$queryPos\t$queryPos\t0\t+\t.\tnote=\"SNP+in+ref
erence:+$refWhat++++Strand+of+query+is+$queryStrand\"\n";
	}
  }

  ## from the coords files we want the regions that are not
  ## covered\n";

  open (F, $fileNameCoords) or die "Problem to open SNP File: $fileNameCoords \n";

  @File=<F>;
  close(F);

  my %coverBB;
  my %coverLB;

  ### get an field, where there is coverage
  foreach (@File) {
	my ($refPos1,$refPos2,$queryPos1,$queryPos2,$overlap1,$overlap2,$identity,$refLength,$queryLength,$dum1,$dum2,$reference,$query) = split(/\s+/);
	$sizeRef{$reference}=$refLength;
	$sizeQuery{$query}=$queryLength;

	$coveredRef{$reference}+=($refPos2-$refPos1+1);

	# check on the reference
	$coverBB{$reference}[$refLength]=undef;
	for ($refPos1..$refPos2) {
	  $coverBB{$reference}[$_]=1;
	}

	# check for the query
	$coverLB{$query}[$queryLength]=undef;
	if ($queryPos2<$queryPos1) {  my $tmp=$queryPos2;$queryPos1=$queryPos2;$queryPos2=$tmp	}
	$coveredQuery{$query}+=(abs($queryPos2-$queryPos1)+1);
	for ($queryPos1..$queryPos2) {
	  $coverLB{$query}[$_]=1;
	}
  }
  ### now parse, were there is no coverage
  for my $chr (keys %coverBB) {
	my $start=1;
	my $ok=1;

	for my $pos (1..(scalar(@{$coverBB{$chr}}))) {
	  if ($ok==1 &&
		  !defined($coverBB{$chr}[$pos])) {
		$ok=0;
		$start=$pos
	  }
	  if ($ok==0 &&
		  defined($coverBB{$chr}[$pos])) {
		$ok=1;
		$BB{$chr} .="$chr\tBBA\tSynteny\t$start\t".($pos-1)."\t0\t+\t.\tnote=\"No s
ynteny with query. Possible insert or too divergent\"\n";
		$noCovRef{$chr}+=($pos-1-$start);

	  }
	}
	if ($ok==0 && $start <scalar(@{$coverBB{$chr}}) ) {
	  $BB{$chr} .="$chr\tBBA\tSynteny\t$start\t".(scalar(@{$coverBB{$chr}})-1)."\t0\t+\
t.\tnote=\"No synteny with query. Possible insert or too divergent\"\n";
	  $noCovRef{$chr}+=($sizeRef{$chr}-1-$start);
	}
  }
  ### now parse, were there is no coverage
  for my $chr (keys %coverLB) {
	my $start=1;
	my $ok=1;

	for my $pos (1..(scalar(@{$coverLB{$chr}})-1)) {
	  if ($ok==1 &&
		  !defined($coverLB{$chr}[$pos])) {
		$ok=0;
		$start=$pos
	  }
	  if ($ok==0 &&
		  defined($coverLB{$chr}[$pos])) {
		$ok=1;
		$LB{$chr} .="$chr\tBBA\tSynteny\t$start\t".($pos-1)."\t0\t+\t.\tnote=\"No s
ynteny with reference. Possible insert or too divergent\"\n";
		$noCovQry{$chr}+=($pos-1-$start);
	  }
	}
	if ($ok==0 && $start < scalar(@{$coverLB{$chr}})) {
	  $LB{$chr} .="$chr\tBBA\tSynteny\t$start\t".(scalar(@{$coverLB{$chr}})-1)."\t0\t+\
t.\tnote=\"No synteny with reference. Possible insert or too divergent\"\n";
	  $noCovQry{$chr}+=($sizeQuery{$chr}-1-$start);
	}
  }
  saveGFF($resultName,"Reference",\%BB);
  saveGFF($resultName,"Query",\%LB);

  my $res;

  foreach my $chr (sort keys %sizeRef ) {
	$res.=sprintf("Of the reference chromosome $chr\t%.2f per cent\thas no synteny with
 the query\n",(( $sizeRef{$chr} - $coveredRef{$chr})*100/$sizeRef{$chr}) );
  }
  foreach my $chr (sort keys %sizeQuery ) {
	$res.=sprintf("Of the query chromosome $chr\t %.2f per cent\thas no synteny with th
e reference\n",((( $sizeQuery{$chr}-$coveredQuery{$chr} ))*100/$sizeQuery{$chr}) );
  }
  return $res
}

############################################
### saveGFF
############################################

sub saveGFF{
  my ($name,$path,$ref_h) = @_;

  if (! -d "$path") {
	!system("mkdir $path") or warn "Couldn't create directory $path.\n";
  }

  foreach my $chr (sort keys %$ref_h){
	my $chrName=$chr;
	$chrName =~ s/\|/_/g;

	open (F,"> $path/$name.$chrName.Mutations.gff") or die "Couldn't create file $name:
 $!\n";

  # UTR must be saved
  print F $$ref_h{$chr};
  close(F);
  }
}
