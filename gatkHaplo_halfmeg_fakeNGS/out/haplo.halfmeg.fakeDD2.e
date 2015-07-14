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
Picked up _JAVA_OPTIONS: -Xmx2000m -Xms50m
INFO  11:55:01,566 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  11:55:01,569 HelpFormatter - The Genome Analysis Toolkit (GATK) v3.4-0-g7e26428, Compiled 2015/05/15 03:25:41 
INFO  11:55:01,569 HelpFormatter - Copyright (c) 2010 The Broad Institute 
INFO  11:55:01,569 HelpFormatter - For support and documentation go to http://www.broadinstitute.org/gatk 
INFO  11:55:01,574 HelpFormatter - Program Args: -T HaplotypeCaller -I bams.DD2.list -ploidy 1 -R /seq/plasmodium/sredmond//refs/PfDD2_v1.fasta -o haplo.rDD2.ALL.vcf -gt_mode DISCOVERY -out_mode EMIT_VARIANTS_ONLY -stand_emit_conf 0 -stand_call_conf 30 -L PfDd2_07_TT:450000-750000 -APO haplo.rDD2.ALL.apo -globalMAPQ -1 -pcrModel NONE -useFilteredReadsForAnnotations -mmq 0 -nda 
INFO  11:55:01,582 HelpFormatter - Executing as sredmond@a2e-bigmem-0.broadinstitute.org on Linux 2.6.32-504.8.1.el6.x86_64 amd64; Java HotSpot(TM) 64-Bit Server VM 1.8.0-b132. 
INFO  11:55:01,582 HelpFormatter - Date/Time: 2015/07/14 11:55:01 
INFO  11:55:01,583 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  11:55:01,583 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  11:55:02,241 GenomeAnalysisEngine - Strictness is SILENT 
INFO  11:55:02,360 GenomeAnalysisEngine - Downsampling Settings: Method: BY_SAMPLE, Target Coverage: 500 
INFO  11:55:02,372 SAMDataSource$SAMReaders - Initializing SAMRecords in serial 
INFO  11:55:02,436 SAMDataSource$SAMReaders - Done initializing BAM readers: total time 0.06 
INFO  11:55:02,458 IntervalUtils - Processing 300001 bp from intervals 
INFO  11:55:02,560 GenomeAnalysisEngine - Preparing for traversal over 3 BAM files 
INFO  11:55:02,601 GenomeAnalysisEngine - Done preparing for traversal 
INFO  11:55:02,601 ProgressMeter - [INITIALIZATION COMPLETE; STARTING PROCESSING] 
INFO  11:55:02,602 ProgressMeter -                 |      processed |    time |         per 1M |           |   total | remaining 
INFO  11:55:02,602 ProgressMeter -        Location | active regions | elapsed | active regions | completed | runtime |   runtime 
INFO  11:55:02,602 HaplotypeCaller - Currently, physical phasing is not available when ploidy is different than 2; therefore it won't be performed 
INFO  11:55:02,725 StrandBiasTest - SAM/BAM data was found. Attempting to use read data to calculate strand bias annotations values. 
INFO  11:55:02,726 StrandBiasTest - SAM/BAM data was found. Attempting to use read data to calculate strand bias annotations values. 
Using SSE4.1 accelerated implementation of PairHMM
INFO  11:55:04,369 VectorLoglessPairHMM - libVectorLoglessPairHMM unpacked successfully from GATK jar file 
INFO  11:55:04,369 VectorLoglessPairHMM - Using vectorized implementation of PairHMM 
WARN  11:55:07,310 HaplotypeScore - Annotation will not be calculated, must be called from UnifiedGenotyper 
WARN  11:55:07,310 InbreedingCoeff - Annotation will not be calculated, must provide a valid PED file (-ped) from the command line. 
INFO  11:55:32,606 ProgressMeter - PfDd2_07_TT:452608              0.0    30.0 s           49.6 w        0.9%    57.5 m      57.0 m 
INFO  11:56:02,607 ProgressMeter - PfDd2_07_TT:454959              0.0    60.0 s           99.2 w        1.7%    60.5 m      59.5 m 
INFO  11:56:32,609 ProgressMeter - PfDd2_07_TT:457215              0.0    90.0 s          148.8 w        2.4%    62.4 m      60.9 m 
WARN  11:56:40,692 AnnotationUtils - Annotation will not be calculated, genotype is not called 
INFO  11:57:02,610 ProgressMeter - PfDd2_07_TT:462510              0.0   120.0 s          198.4 w        4.2%    48.0 m      46.0 m 
INFO  11:57:32,612 ProgressMeter - PfDd2_07_TT:467590              0.0     2.5 m          248.0 w        5.9%    42.6 m      40.1 m 
INFO  11:58:02,613 ProgressMeter - PfDd2_07_TT:469884              0.0     3.0 m          297.6 w        6.6%    45.3 m      42.3 m 
INFO  11:58:32,615 ProgressMeter - PfDd2_07_TT:472109              0.0     3.5 m          347.2 w        7.4%    47.5 m      44.0 m 
INFO  11:59:02,617 ProgressMeter - PfDd2_07_TT:474456              0.0     4.0 m          396.9 w        8.2%    49.1 m      45.1 m 
INFO  11:59:32,618 ProgressMeter - PfDd2_07_TT:478911              0.0     4.5 m          446.5 w        9.6%    46.7 m      42.2 m 
INFO  12:00:02,620 ProgressMeter - PfDd2_07_TT:480252              0.0     5.0 m          496.1 w       10.1%    49.6 m      44.6 m 
INFO  12:00:32,622 ProgressMeter - PfDd2_07_TT:483545              0.0     5.5 m          545.7 w       11.2%    49.2 m      43.7 m 
INFO  12:01:02,623 ProgressMeter - PfDd2_07_TT:485199              0.0     6.0 m          595.3 w       11.7%    51.1 m      45.1 m 
INFO  12:01:32,625 ProgressMeter - PfDd2_07_TT:486430              0.0     6.5 m          644.9 w       12.1%    53.5 m      47.0 m 
INFO  12:02:02,626 ProgressMeter - PfDd2_07_TT:488353              0.0     7.0 m          694.5 w       12.8%    54.8 m      47.8 m 
INFO  12:02:32,628 ProgressMeter - PfDd2_07_TT:496636              0.0     7.5 m          744.1 w       15.5%    48.2 m      40.7 m 
INFO  12:03:02,631 ProgressMeter - PfDd2_07_TT:503280              0.0     8.0 m          793.7 w       17.8%    45.0 m      37.0 m 
INFO  12:03:32,633 ProgressMeter - PfDd2_07_TT:509730              0.0     8.5 m          843.3 w       19.9%    42.7 m      34.2 m 
INFO  12:04:02,634 ProgressMeter - PfDd2_07_TT:516686              0.0     9.0 m          892.9 w       22.2%    40.5 m      31.5 m 
INFO  12:04:32,636 ProgressMeter - PfDd2_07_TT:518853              0.0     9.5 m          942.5 w       23.0%    41.4 m      31.9 m 
INFO  12:05:02,637 ProgressMeter - PfDd2_07_TT:521824              0.0    10.0 m          992.1 w       23.9%    41.8 m      31.8 m 
INFO  12:05:32,639 ProgressMeter - PfDd2_07_TT:528883              0.0    10.5 m         1041.7 w       26.3%    39.9 m      29.4 m 
INFO  12:06:02,640 ProgressMeter - PfDd2_07_TT:536598              0.0    11.0 m         1091.3 w       28.9%    38.1 m      27.1 m 
INFO  12:06:32,642 ProgressMeter - PfDd2_07_TT:550363              0.0    11.5 m         1140.9 w       33.5%    34.4 m      22.9 m 
INFO  12:07:02,643 ProgressMeter - PfDd2_07_TT:552645              0.0    12.0 m         1190.5 w       34.2%    35.1 m      23.1 m 
INFO  12:07:32,645 ProgressMeter - PfDd2_07_TT:554908              0.0    12.5 m         1240.2 w       35.0%    35.7 m      23.2 m 
INFO  12:08:02,646 ProgressMeter - PfDd2_07_TT:557089              0.0    13.0 m         1289.8 w       35.7%    36.4 m      23.4 m 
INFO  12:08:32,648 ProgressMeter - PfDd2_07_TT:559919              0.0    13.5 m         1339.4 w       36.6%    36.8 m      23.3 m 
INFO  12:09:02,649 ProgressMeter - PfDd2_07_TT:561047              0.0    14.0 m         1389.0 w       37.0%    37.8 m      23.8 m 
INFO  12:09:32,651 ProgressMeter - PfDd2_07_TT:562234              0.0    14.5 m         1438.6 w       37.4%    38.8 m      24.3 m 
INFO  12:10:02,652 ProgressMeter - PfDd2_07_TT:563208              0.0    15.0 m         1488.2 w       37.7%    39.7 m      24.7 m 
INFO  12:10:32,653 ProgressMeter - PfDd2_07_TT:564741              0.0    15.5 m         1537.8 w       38.2%    40.5 m      25.0 m 
INFO  12:11:02,655 ProgressMeter - PfDd2_07_TT:565452              0.0    16.0 m         1587.4 w       38.5%    41.6 m      25.6 m 
INFO  12:11:32,656 ProgressMeter - PfDd2_07_TT:566479              0.0    16.5 m         1637.0 w       38.8%    42.5 m      26.0 m 
INFO  12:12:02,658 ProgressMeter - PfDd2_07_TT:569230              0.0    17.0 m         1686.6 w       39.7%    42.8 m      25.8 m 
INFO  12:12:32,659 ProgressMeter - PfDd2_07_TT:570843              0.0    17.5 m         1736.2 w       40.3%    43.4 m      25.9 m 
INFO  12:13:02,661 ProgressMeter - PfDd2_07_TT:573046              0.0    18.0 m         1785.8 w       41.0%    43.9 m      25.9 m 
INFO  12:13:32,671 ProgressMeter - PfDd2_07_TT:575652              0.0    18.5 m         1835.4 w       41.9%    44.2 m      25.7 m 
INFO  12:14:02,672 ProgressMeter - PfDd2_07_TT:578663              0.0    19.0 m         1885.0 w       42.9%    44.3 m      25.3 m 
INFO  12:14:32,674 ProgressMeter - PfDd2_07_TT:582849              0.0    19.5 m         1934.6 w       44.3%    44.0 m      24.5 m 
INFO  12:15:02,675 ProgressMeter - PfDd2_07_TT:587720              0.0    20.0 m         1984.2 w       45.9%    43.6 m      23.6 m 
INFO  12:15:32,677 ProgressMeter - PfDd2_07_TT:590343              0.0    20.5 m         2033.9 w       46.8%    43.8 m      23.3 m 
INFO  12:16:02,679 ProgressMeter - PfDd2_07_TT:595515              0.0    21.0 m         2083.5 w       48.5%    43.3 m      22.3 m 
INFO  12:16:32,681 ProgressMeter - PfDd2_07_TT:599685              0.0    21.5 m         2133.1 w       49.9%    43.1 m      21.6 m 
INFO  12:17:02,682 ProgressMeter - PfDd2_07_TT:603725              0.0    22.0 m         2182.7 w       51.2%    42.9 m      20.9 m 
INFO  12:17:32,684 ProgressMeter - PfDd2_07_TT:606757              0.0    22.5 m         2232.3 w       52.3%    43.1 m      20.6 m 
INFO  12:18:02,685 ProgressMeter - PfDd2_07_TT:611159              0.0    23.0 m         2281.9 w       53.7%    42.8 m      19.8 m 
INFO  12:18:32,687 ProgressMeter - PfDd2_07_TT:616115              0.0    23.5 m         2331.5 w       55.4%    42.4 m      18.9 m 
INFO  12:19:02,688 ProgressMeter - PfDd2_07_TT:624521              0.0    24.0 m         2381.1 w       58.2%    41.3 m      17.3 m 
INFO  12:19:32,690 ProgressMeter - PfDd2_07_TT:627435              0.0    24.5 m         2430.7 w       59.1%    41.4 m      16.9 m 
INFO  12:20:02,692 ProgressMeter - PfDd2_07_TT:631033              0.0    25.0 m         2480.3 w       60.3%    41.4 m      16.4 m 
INFO  12:20:32,693 ProgressMeter - PfDd2_07_TT:634028              0.0    25.5 m         2529.9 w       61.3%    41.6 m      16.1 m 
INFO  12:21:02,694 ProgressMeter - PfDd2_07_TT:636993              0.0    26.0 m         2579.5 w       62.3%    41.7 m      15.7 m 
INFO  12:21:32,696 ProgressMeter - PfDd2_07_TT:642191              0.0    26.5 m         2629.1 w       64.1%    41.4 m      14.9 m 
INFO  12:22:02,697 ProgressMeter - PfDd2_07_TT:645741              0.0    27.0 m         2678.7 w       65.2%    41.4 m      14.4 m 
INFO  12:22:32,699 ProgressMeter - PfDd2_07_TT:649042              0.0    27.5 m         2728.3 w       66.3%    41.4 m      13.9 m 
INFO  12:23:02,701 ProgressMeter - PfDd2_07_TT:651892              0.0    28.0 m         2777.9 w       67.3%    41.6 m      13.6 m 
INFO  12:23:32,702 ProgressMeter - PfDd2_07_TT:656147              0.0    28.5 m         2827.5 w       68.7%    41.5 m      13.0 m 
INFO  12:24:02,726 ProgressMeter - PfDd2_07_TT:658586              0.0    29.0 m         2877.2 w       69.5%    41.7 m      12.7 m 
INFO  12:24:32,727 ProgressMeter - PfDd2_07_TT:662153              0.0    29.5 m         2926.8 w       70.7%    41.7 m      12.2 m 
INFO  12:25:02,729 ProgressMeter - PfDd2_07_TT:665553              0.0    30.0 m         2976.4 w       71.9%    41.8 m      11.8 m 
INFO  12:25:32,730 ProgressMeter - PfDd2_07_TT:668214              0.0    30.5 m         3026.0 w       72.7%    41.9 m      11.4 m 
WARN  12:26:01,222 ExactAFCalculator - this tool is currently set to genotype at most 6 alternate alleles in a given context, but the context at PfDd2_07_TT:668985 has 11 alternate alleles so only the top alleles will be used; see the --max_alternate_alleles argument 
INFO  12:26:02,731 ProgressMeter - PfDd2_07_TT:669467              0.0    31.0 m         3075.6 w       73.2%    42.4 m      11.4 m 
INFO  12:26:32,733 ProgressMeter - PfDd2_07_TT:674005              0.0    31.5 m         3125.2 w       74.7%    42.2 m      10.7 m 
INFO  12:27:02,734 ProgressMeter - PfDd2_07_TT:677014              0.0    32.0 m         3174.8 w       75.7%    42.3 m      10.3 m 
INFO  12:27:32,748 ProgressMeter - PfDd2_07_TT:678548              0.0    32.5 m         3224.4 w       76.2%    42.7 m      10.2 m 
WARN  12:27:38,901 ExactAFCalculator - this tool is currently set to genotype at most 6 alternate alleles in a given context, but the context at PfDd2_07_TT:678235 has 10 alternate alleles so only the top alleles will be used; see the --max_alternate_alleles argument 
INFO  12:28:02,749 ProgressMeter - PfDd2_07_TT:680890              0.0    33.0 m         3274.1 w       77.0%    42.9 m       9.9 m 
INFO  12:28:32,750 ProgressMeter - PfDd2_07_TT:683305              0.0    33.5 m         3323.7 w       77.8%    43.1 m       9.6 m 
INFO  12:29:02,752 ProgressMeter - PfDd2_07_TT:685527              0.0    34.0 m         3373.3 w       78.5%    43.3 m       9.3 m 
WARN  12:29:06,905 ExactAFCalculator - this tool is currently set to genotype at most 6 alternate alleles in a given context, but the context at PfDd2_07_TT:685223 has 8 alternate alleles so only the top alleles will be used; see the --max_alternate_alleles argument 
INFO  12:29:32,753 ProgressMeter - PfDd2_07_TT:687876              0.0    34.5 m         3422.9 w       79.3%    43.5 m       9.0 m 
INFO  12:30:02,755 ProgressMeter - PfDd2_07_TT:689746              0.0    35.0 m         3472.5 w       79.9%    43.8 m       8.8 m 
INFO  12:30:32,756 ProgressMeter - PfDd2_07_TT:691951              0.0    35.5 m         3522.1 w       80.7%    44.0 m       8.5 m 
INFO  12:31:02,758 ProgressMeter - PfDd2_07_TT:693913              0.0    36.0 m         3571.7 w       81.3%    44.3 m       8.3 m 
INFO  12:31:32,761 ProgressMeter - PfDd2_07_TT:695768              0.0    36.5 m         3621.3 w       81.9%    44.6 m       8.1 m 
INFO  12:32:02,762 ProgressMeter - PfDd2_07_TT:698359              0.0    37.0 m         3670.9 w       82.8%    44.7 m       7.7 m 
INFO  12:32:32,764 ProgressMeter - PfDd2_07_TT:700296              0.0    37.5 m         3720.5 w       83.4%    44.9 m       7.4 m 
INFO  12:33:02,765 ProgressMeter - PfDd2_07_TT:702028              0.0    38.0 m         3770.1 w       84.0%    45.2 m       7.2 m 
INFO  12:33:32,767 ProgressMeter - PfDd2_07_TT:705483              0.0    38.5 m         3819.7 w       85.2%    45.2 m       6.7 m 
INFO  12:34:02,769 ProgressMeter - PfDd2_07_TT:707935              0.0    39.0 m         3869.3 w       86.0%    45.4 m       6.4 m 
INFO  12:34:32,770 ProgressMeter - PfDd2_07_TT:711380              0.0    39.5 m         3918.9 w       87.1%    45.3 m       5.8 m 
INFO  12:35:02,784 ProgressMeter - PfDd2_07_TT:717132              0.0    40.0 m         3968.6 w       89.0%    44.9 m       4.9 m 
INFO  12:35:32,785 ProgressMeter - PfDd2_07_TT:723630              0.0    40.5 m         4018.2 w       91.2%    44.4 m       3.9 m 
INFO  12:35:43,678 VectorLoglessPairHMM - Time spent in setup for JNI call : 0.33510401100000003 
INFO  12:35:43,679 PairHMM - Total compute time in PairHMM computeLikelihoods() : 1977.675293938 
INFO  12:35:43,680 HaplotypeCaller - Ran local assembly on 566 active regions 
INFO  12:35:43,706 ProgressMeter -            done         300001.0    40.7 m            2.3 h      100.0%    40.7 m       0.0 s 
INFO  12:35:43,707 ProgressMeter - Total runtime 2441.11 secs, 40.69 min, 0.68 hours 
INFO  12:35:43,708 MicroScheduler - 2274 reads were filtered out during the traversal out of approximately 293730 total reads (0.77%) 
INFO  12:35:43,708 MicroScheduler -   -> 0 reads (0.00% of total) failing BadCigarFilter 
INFO  12:35:43,708 MicroScheduler -   -> 0 reads (0.00% of total) failing DuplicateReadFilter 
INFO  12:35:43,709 MicroScheduler -   -> 0 reads (0.00% of total) failing FailsVendorQualityCheckFilter 
INFO  12:35:43,709 MicroScheduler -   -> 0 reads (0.00% of total) failing HCMappingQualityFilter 
INFO  12:35:43,709 MicroScheduler -   -> 0 reads (0.00% of total) failing MalformedReadFilter 
INFO  12:35:43,712 MicroScheduler -   -> 0 reads (0.00% of total) failing MappingQualityUnavailableFilter 
INFO  12:35:43,712 MicroScheduler -   -> 2102 reads (0.72% of total) failing NotPrimaryAlignmentFilter 
INFO  12:35:43,712 MicroScheduler -   -> 172 reads (0.06% of total) failing UnmappedReadFilter 
INFO  12:35:44,549 GATKRunReport - Uploaded run statistics report to AWS S3 
