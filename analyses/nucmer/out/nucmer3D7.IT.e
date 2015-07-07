1: PREPARING DATA
2,3: RUNNING mummer AND CREATING CLUSTERS
# reading input file "nucmer.3D7vIT.ntref" of length 23332847
# construct suffix tree for sequence of length 23332847
# (maximum reference length is 536870908)
# (maximum query length is 4294967295)
# process 233328 characters per dot
#....................................................................................................
# CONSTRUCTIONTIME /broad/software/free/Linux/redhat_6_x86_64/pkgs/mummer_3.22/mummer nucmer.3D7vIT.ntref 14.70
# reading input file "/seq/plasmodium/sredmond/pfdisco/analyses/nucmer/PfIT_v3.fasta" of length 22983634
# matching query-file "/seq/plasmodium/sredmond/pfdisco/analyses/nucmer/PfIT_v3.fasta"
# against subject-file "nucmer.3D7vIT.ntref"
# COMPLETETIME /broad/software/free/Linux/redhat_6_x86_64/pkgs/mummer_3.22/mummer nucmer.3D7vIT.ntref 44.28
# SPACE /broad/software/free/Linux/redhat_6_x86_64/pkgs/mummer_3.22/mummer nucmer.3D7vIT.ntref 44.91
4: FINISHING DATA
WARNING: Invalid break length 20000, ignoring
5: GENERATING COORDS FILE
/seq/plasmodium/sredmond/pfdisco/analyses/nucmer/Pf3D7_v3.fasta	/seq/plasmodium/sredmond/pfdisco/analyses/nucmer/PfIT_v3.fasta
nucmer.Pf3D7_v3xPfIT_v3.bed	nucmer.PfIT_v3xPf3D7_v3.bed
Traceback (most recent call last):
  File "/seq/plasmodium/sredmond//pfdisco//scripts/mummer2vcf_plus1error.py", line 359, in <module>
    _write_var()
  File "/seq/plasmodium/sredmond//pfdisco//scripts/mummer2vcf_plus1error.py", line 104, in _write_var
    seq2_2 = pb2+revcomp(seq2)
  File "/seq/plasmodium/sredmond//pfdisco//scripts/mummer2vcf_plus1error.py", line 219, in revcomp
    r = [d[b] for b in reversed(seq)]
KeyError: 'N'
