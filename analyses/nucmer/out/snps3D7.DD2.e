  File "/seq/plasmodium/sredmond//pfdisco//scripts/mummer2vcf.py", line 316
    vcfdict[sample1]=vcf1
                        ^
IndentationError: unindent does not match any outer indentation level
/seq/plasmodium/sredmond//pfdisco//scripts/runNucmerPipeline.sh: line 8: nucmer: command not found
/seq/plasmodium/sredmond//pfdisco//scripts/runNucmerPipeline.sh: line 11: show-snps: command not found
Traceback (most recent call last):
  File "/seq/plasmodium/sredmond//pfdisco//scripts/mummer2vcf_plus1error.py", line 371, in <module>
    _write_var()
  File "/seq/plasmodium/sredmond//pfdisco//scripts/mummer2vcf_plus1error.py", line 76, in _write_var
    pb1 = _get_base(fasta1,lc1,st1-1)
  File "/seq/plasmodium/sredmond//pfdisco//scripts/mummer2vcf_plus1error.py", line 213, in _get_base
    prevbase = getseq(fasta, chrom, pos)
  File "/seq/plasmodium/sredmond//pfdisco//scripts/mummer2vcf_plus1error.py", line 224, in getseq
    command = ['samtools', 'faidx', fasta, chrom+":"+str(pos)+"-"+str(pos2)]
TypeError: unsupported operand type(s) for +: 'NoneType' and 'str'
/seq/plasmodium/sredmond//pfdisco//scripts/runNucmerPipeline.sh: line 8: nucmer: command not found
/seq/plasmodium/sredmond//pfdisco//scripts/runNucmerPipeline.sh: line 11: show-snps: command not found
Traceback (most recent call last):
  File "/seq/plasmodium/sredmond//pfdisco//scripts/mummer2vcf_plus1error.py", line 371, in <module>
    _write_var()
  File "/seq/plasmodium/sredmond//pfdisco//scripts/mummer2vcf_plus1error.py", line 76, in _write_var
    pb1 = _get_base(fasta1,lc1,st1-1)
  File "/seq/plasmodium/sredmond//pfdisco//scripts/mummer2vcf_plus1error.py", line 213, in _get_base
    prevbase = getseq(fasta, chrom, pos)
  File "/seq/plasmodium/sredmond//pfdisco//scripts/mummer2vcf_plus1error.py", line 224, in getseq
    command = ['samtools', 'faidx', fasta, chrom+":"+str(pos)+"-"+str(pos2)]
TypeError: unsupported operand type(s) for +: 'NoneType' and 'str'
/seq/plasmodium/sredmond//pfdisco//scripts/runNucmerPipeline.sh: line 8: nucmer: command not found
/seq/plasmodium/sredmond//pfdisco//scripts/runNucmerPipeline.sh: line 11: show-snps: command not found
Traceback (most recent call last):
  File "/seq/plasmodium/sredmond//pfdisco//scripts/mummer2vcf_plus1error.py", line 371, in <module>
    _write_var()
  File "/seq/plasmodium/sredmond//pfdisco//scripts/mummer2vcf_plus1error.py", line 76, in _write_var
    pb1 = _get_base(fasta1,lc1,st1-1)
  File "/seq/plasmodium/sredmond//pfdisco//scripts/mummer2vcf_plus1error.py", line 213, in _get_base
    prevbase = getseq(fasta, chrom, pos)
  File "/seq/plasmodium/sredmond//pfdisco//scripts/mummer2vcf_plus1error.py", line 224, in getseq
    command = ['samtools', 'faidx', fasta, chrom+":"+str(pos)+"-"+str(pos2)]
TypeError: unsupported operand type(s) for +: 'NoneType' and 'str'
1: PREPARING DATA
2,3: RUNNING mummer AND CREATING CLUSTERS
# reading input file "nucmer.3D7vIT.ntref" of length 23332847
# construct suffix tree for sequence of length 23332847
# (maximum reference length is 536870908)
# (maximum query length is 4294967295)
# process 233328 characters per dot
#....................................................................................................
# CONSTRUCTIONTIME /broad/software/free/Linux/redhat_6_x86_64/pkgs/mummer_3.22/mummer nucmer.3D7vIT.ntref 14.27
# reading input file "/seq/plasmodium/sredmond/pfdisco/analyses/nucmer/PfIT_v3.fasta" of length 22983634
# matching query-file "/seq/plasmodium/sredmond/pfdisco/analyses/nucmer/PfIT_v3.fasta"
# against subject-file "nucmer.3D7vIT.ntref"
# COMPLETETIME /broad/software/free/Linux/redhat_6_x86_64/pkgs/mummer_3.22/mummer nucmer.3D7vIT.ntref 46.10
# SPACE /broad/software/free/Linux/redhat_6_x86_64/pkgs/mummer_3.22/mummer nucmer.3D7vIT.ntref 44.91
4: FINISHING DATA
WARNING: Invalid break length 20000, ignoring
5: GENERATING COORDS FILE
/seq/plasmodium/sredmond/pfdisco/analyses/nucmer/Pf3D7_v3.fasta	/seq/plasmodium/sredmond/pfdisco/analyses/nucmer/PfIT_v3.fasta
nucmer.Pf3D7_v3inPfIT_v3.bed	nucmer.PfIT_v3inPf3D7_v3.bed
Traceback (most recent call last):
  File "/seq/plasmodium/sredmond//pfdisco//scripts/mummer2vcf_plus1error.py", line 359, in <module>
    _write_var()
  File "/seq/plasmodium/sredmond//pfdisco//scripts/mummer2vcf_plus1error.py", line 104, in _write_var
    seq2_2 = pb2+revcomp(seq2)
  File "/seq/plasmodium/sredmond//pfdisco//scripts/mummer2vcf_plus1error.py", line 219, in revcomp
    r = [d[b] for b in reversed(seq)]
KeyError: 'N'
