Picked up _JAVA_OPTIONS: -Xmx2000m -Xms50m
INFO  20:00:03,774 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  20:00:03,777 HelpFormatter - The Genome Analysis Toolkit (GATK) v3.4-0-g7e26428, Compiled 2015/05/15 03:25:41 
INFO  20:00:03,778 HelpFormatter - Copyright (c) 2010 The Broad Institute 
INFO  20:00:03,778 HelpFormatter - For support and documentation go to http://www.broadinstitute.org/gatk 
INFO  20:00:03,783 HelpFormatter - Program Args: -T HaplotypeCaller -I bams.DD2.list -ploidy 1 -R /seq/plasmodium/sredmond//refs/PfDD2_v1.fasta -o haplo.rDD2.ALL.vcf -gt_mode DISCOVERY -out_mode EMIT_VARIANTS_ONLY -stand_emit_conf 0 -stand_call_conf 30 -L PfDD2_07_v1:450000-750000 -APO haplo.rDD2.ALL.apo -globalMAPQ -1 -pcrModel NONE -useFilteredReadsForAnnotations -mmq 0 -nda 
INFO  20:00:03,791 HelpFormatter - Executing as sredmond@a2e-bigmem-0.broadinstitute.org on Linux 2.6.32-504.8.1.el6.x86_64 amd64; Java HotSpot(TM) 64-Bit Server VM 1.8.0-b132. 
INFO  20:00:03,791 HelpFormatter - Date/Time: 2015/07/13 20:00:03 
INFO  20:00:03,792 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  20:00:03,792 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  20:00:04,479 GenomeAnalysisEngine - Strictness is SILENT 
INFO  20:00:04,768 GenomeAnalysisEngine - Downsampling Settings: Method: BY_SAMPLE, Target Coverage: 500 
INFO  20:00:04,785 SAMDataSource$SAMReaders - Initializing SAMRecords in serial 
INFO  20:00:05,532 SAMDataSource$SAMReaders - Done initializing BAM readers: total time 0.75 
INFO  20:00:06,318 GATKRunReport - Uploaded run statistics report to AWS S3 
##### ERROR ------------------------------------------------------------------------------------------
##### ERROR A USER ERROR has occurred (version 3.4-0-g7e26428): 
##### ERROR
##### ERROR This means that one or more arguments or inputs in your command are incorrect.
##### ERROR The error message below tells you what is the problem.
##### ERROR
##### ERROR If the problem is an invalid argument, please check the online documentation guide
##### ERROR (or rerun your command with --help) to view allowable command-line arguments for this tool.
##### ERROR
##### ERROR Visit our website and forum for extensive documentation and answers to 
##### ERROR commonly asked questions http://www.broadinstitute.org/gatk
##### ERROR
##### ERROR Please do NOT post this error to the GATK forum unless you have really tried to fix it yourself.
##### ERROR
##### ERROR MESSAGE: Badly formed genome loc: Contig 'PfDD2_07_v1' does not match any contig in the GATK sequence dictionary derived from the reference; are you sure you are using the correct reference fasta file?
##### ERROR ------------------------------------------------------------------------------------------
Picked up _JAVA_OPTIONS: -Xmx2000m -Xms50m
INFO  20:40:46,287 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  20:40:46,291 HelpFormatter - The Genome Analysis Toolkit (GATK) v3.4-0-g7e26428, Compiled 2015/05/15 03:25:41 
INFO  20:40:46,292 HelpFormatter - Copyright (c) 2010 The Broad Institute 
INFO  20:40:46,292 HelpFormatter - For support and documentation go to http://www.broadinstitute.org/gatk 
INFO  20:40:46,299 HelpFormatter - Program Args: -T HaplotypeCaller -I bams.DD2.list -ploidy 1 -R /seq/plasmodium/sredmond//refs/PfDD2_v1.fasta -o haplo.rDD2.ALL.vcf -gt_mode DISCOVERY -out_mode EMIT_VARIANTS_ONLY -stand_emit_conf 0 -stand_call_conf 30 -L PfDD2_07_TT:450000-750000 -APO haplo.rDD2.ALL.apo -globalMAPQ -1 -pcrModel NONE -useFilteredReadsForAnnotations -mmq 0 -nda 
INFO  20:40:46,308 HelpFormatter - Executing as sredmond@a2e-bigmem-0.broadinstitute.org on Linux 2.6.32-504.8.1.el6.x86_64 amd64; Java HotSpot(TM) 64-Bit Server VM 1.8.0-b132. 
INFO  20:40:46,308 HelpFormatter - Date/Time: 2015/07/13 20:40:46 
INFO  20:40:46,308 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  20:40:46,309 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  20:40:47,062 GenomeAnalysisEngine - Strictness is SILENT 
INFO  20:40:47,333 GenomeAnalysisEngine - Downsampling Settings: Method: BY_SAMPLE, Target Coverage: 500 
INFO  20:40:47,343 SAMDataSource$SAMReaders - Initializing SAMRecords in serial 
INFO  20:40:47,408 SAMDataSource$SAMReaders - Done initializing BAM readers: total time 0.06 
INFO  20:40:48,258 GATKRunReport - Uploaded run statistics report to AWS S3 
##### ERROR ------------------------------------------------------------------------------------------
##### ERROR A USER ERROR has occurred (version 3.4-0-g7e26428): 
##### ERROR
##### ERROR This means that one or more arguments or inputs in your command are incorrect.
##### ERROR The error message below tells you what is the problem.
##### ERROR
##### ERROR If the problem is an invalid argument, please check the online documentation guide
##### ERROR (or rerun your command with --help) to view allowable command-line arguments for this tool.
##### ERROR
##### ERROR Visit our website and forum for extensive documentation and answers to 
##### ERROR commonly asked questions http://www.broadinstitute.org/gatk
##### ERROR
##### ERROR Please do NOT post this error to the GATK forum unless you have really tried to fix it yourself.
##### ERROR
##### ERROR MESSAGE: Badly formed genome loc: Contig 'PfDD2_07_TT' does not match any contig in the GATK sequence dictionary derived from the reference; are you sure you are using the correct reference fasta file?
##### ERROR ------------------------------------------------------------------------------------------
Picked up _JAVA_OPTIONS: -Xmx2000m -Xms50m
INFO  11:53:24,727 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  11:53:24,731 HelpFormatter - The Genome Analysis Toolkit (GATK) v3.4-0-g7e26428, Compiled 2015/05/15 03:25:41 
INFO  11:53:24,731 HelpFormatter - Copyright (c) 2010 The Broad Institute 
INFO  11:53:24,732 HelpFormatter - For support and documentation go to http://www.broadinstitute.org/gatk 
INFO  11:53:24,740 HelpFormatter - Program Args: -T HaplotypeCaller -I bams.DD2.list -ploidy 1 -R /seq/plasmodium/sredmond//refs/PfDD2_v1.fasta -o haplo.rDD2.ALL.vcf -gt_mode DISCOVERY -out_mode EMIT_VARIANTS_ONLY -stand_emit_conf 0 -stand_call_conf 30 -L PfDD2_07_TT:450000-750000 -APO haplo.rDD2.ALL.apo -globalMAPQ -1 -pcrModel NONE -useFilteredReadsForAnnotations -mmq 0 -nda 
INFO  11:53:24,748 HelpFormatter - Executing as sredmond@a2e-bigmem-0.broadinstitute.org on Linux 2.6.32-504.8.1.el6.x86_64 amd64; Java HotSpot(TM) 64-Bit Server VM 1.8.0-b132. 
INFO  11:53:24,748 HelpFormatter - Date/Time: 2015/07/14 11:53:24 
INFO  11:53:24,748 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  11:53:24,748 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  11:53:25,394 GenomeAnalysisEngine - Strictness is SILENT 
INFO  11:53:25,668 GenomeAnalysisEngine - Downsampling Settings: Method: BY_SAMPLE, Target Coverage: 500 
INFO  11:53:25,678 SAMDataSource$SAMReaders - Initializing SAMRecords in serial 
INFO  11:53:25,890 SAMDataSource$SAMReaders - Done initializing BAM readers: total time 0.21 
INFO  11:53:26,614 GATKRunReport - Uploaded run statistics report to AWS S3 
##### ERROR ------------------------------------------------------------------------------------------
##### ERROR A USER ERROR has occurred (version 3.4-0-g7e26428): 
##### ERROR
##### ERROR This means that one or more arguments or inputs in your command are incorrect.
##### ERROR The error message below tells you what is the problem.
##### ERROR
##### ERROR If the problem is an invalid argument, please check the online documentation guide
##### ERROR (or rerun your command with --help) to view allowable command-line arguments for this tool.
##### ERROR
##### ERROR Visit our website and forum for extensive documentation and answers to 
##### ERROR commonly asked questions http://www.broadinstitute.org/gatk
##### ERROR
##### ERROR Please do NOT post this error to the GATK forum unless you have really tried to fix it yourself.
##### ERROR
##### ERROR MESSAGE: Badly formed genome loc: Contig 'PfDD2_07_TT' does not match any contig in the GATK sequence dictionary derived from the reference; are you sure you are using the correct reference fasta file?
##### ERROR ------------------------------------------------------------------------------------------
