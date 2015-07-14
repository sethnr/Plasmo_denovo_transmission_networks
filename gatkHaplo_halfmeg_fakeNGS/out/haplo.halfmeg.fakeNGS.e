Picked up _JAVA_OPTIONS: -Xmx2000m -Xms50m
INFO  13:52:42,566 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  13:52:42,569 HelpFormatter - The Genome Analysis Toolkit (GATK) v3.4-0-g7e26428, Compiled 2015/05/15 03:25:41 
INFO  13:52:42,569 HelpFormatter - Copyright (c) 2010 The Broad Institute 
INFO  13:52:42,569 HelpFormatter - For support and documentation go to http://www.broadinstitute.org/gatk 
INFO  13:52:42,574 HelpFormatter - Program Args: -T HaplotypeCaller -I bams.3D7.list -ploidy 1 -R /seq/plasmodium/sredmond//refs/Pf3D7_v3.fasta -o haplo.r3D7.ALL.vcf -gt_mode DISCOVERY -out_mode EMIT_VARIANTS_ONLY -stand_emit_conf 0 -stand_call_conf 30 -L Pf3D7_07_v3:450000-750000 -APO haplo.r3D7.ALL.apo -globalMAPQ -1 -pcrModel NONE -useFilteredReadsForAnnotations -mmq 0 -nda 
INFO  13:52:42,582 HelpFormatter - Executing as sredmond@a2e-bigmem-0.broadinstitute.org on Linux 2.6.32-504.8.1.el6.x86_64 amd64; Java HotSpot(TM) 64-Bit Server VM 1.8.0-b132. 
INFO  13:52:42,582 HelpFormatter - Date/Time: 2015/07/13 13:52:42 
INFO  13:52:42,582 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  13:52:42,583 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  13:52:43,449 GenomeAnalysisEngine - Strictness is SILENT 
INFO  13:52:43,566 GenomeAnalysisEngine - Downsampling Settings: Method: BY_SAMPLE, Target Coverage: 500 
INFO  13:52:43,577 SAMDataSource$SAMReaders - Initializing SAMRecords in serial 
INFO  13:52:43,768 SAMDataSource$SAMReaders - Done initializing BAM readers: total time 0.19 
INFO  13:52:43,811 IntervalUtils - Processing 300001 bp from intervals 
INFO  13:52:43,910 GenomeAnalysisEngine - Preparing for traversal over 3 BAM files 
INFO  13:52:43,945 GenomeAnalysisEngine - Done preparing for traversal 
INFO  13:52:43,946 ProgressMeter - [INITIALIZATION COMPLETE; STARTING PROCESSING] 
INFO  13:52:43,947 ProgressMeter -                 |      processed |    time |         per 1M |           |   total | remaining 
INFO  13:52:43,947 ProgressMeter -        Location | active regions | elapsed | active regions | completed | runtime |   runtime 
INFO  13:52:43,947 HaplotypeCaller - Currently, physical phasing is not available when ploidy is different than 2; therefore it won't be performed 
INFO  13:52:44,049 StrandBiasTest - SAM/BAM data was found. Attempting to use read data to calculate strand bias annotations values. 
INFO  13:52:44,049 StrandBiasTest - SAM/BAM data was found. Attempting to use read data to calculate strand bias annotations values. 
Using SSE4.1 accelerated implementation of PairHMM
INFO  13:52:45,690 VectorLoglessPairHMM - libVectorLoglessPairHMM unpacked successfully from GATK jar file 
INFO  13:52:45,691 VectorLoglessPairHMM - Using vectorized implementation of PairHMM 
Picked up _JAVA_OPTIONS: -Xmx2000m -Xms50m
INFO  13:53:18,449 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  13:53:18,452 HelpFormatter - The Genome Analysis Toolkit (GATK) v3.4-0-g7e26428, Compiled 2015/05/15 03:25:41 
INFO  13:53:18,452 HelpFormatter - Copyright (c) 2010 The Broad Institute 
INFO  13:53:18,452 HelpFormatter - For support and documentation go to http://www.broadinstitute.org/gatk 
INFO  13:53:18,457 HelpFormatter - Program Args: -T HaplotypeCaller -I bams.3D7.list -ploidy 1 -R /seq/plasmodium/sredmond//refs/Pf3D7_v3.fasta -o haplo.r3D7.ALL.48-50.vcf -gt_mode DISCOVERY -out_mode EMIT_VARIANTS_ONLY -stand_emit_conf 0 -stand_call_conf 30 -L Pf3D7_07_v3:480000-500000 -APO haplo.r3D7.ALL.48-50.apo -globalMAPQ -1 -pcrModel NONE -useFilteredReadsForAnnotations -mmq 0 -nda 
INFO  13:53:18,464 HelpFormatter - Executing as sredmond@a2e-bigmem-0.broadinstitute.org on Linux 2.6.32-504.8.1.el6.x86_64 amd64; Java HotSpot(TM) 64-Bit Server VM 1.8.0-b132. 
INFO  13:53:18,464 HelpFormatter - Date/Time: 2015/07/13 13:53:18 
INFO  13:53:18,464 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  13:53:18,464 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  13:53:19,122 GenomeAnalysisEngine - Strictness is SILENT 
INFO  13:53:19,237 GenomeAnalysisEngine - Downsampling Settings: Method: BY_SAMPLE, Target Coverage: 500 
INFO  13:53:19,250 SAMDataSource$SAMReaders - Initializing SAMRecords in serial 
INFO  13:53:19,309 SAMDataSource$SAMReaders - Done initializing BAM readers: total time 0.06 
INFO  13:53:19,330 IntervalUtils - Processing 20001 bp from intervals 
INFO  13:53:19,438 GenomeAnalysisEngine - Preparing for traversal over 3 BAM files 
INFO  13:53:19,548 GenomeAnalysisEngine - Done preparing for traversal 
INFO  13:53:19,550 ProgressMeter - [INITIALIZATION COMPLETE; STARTING PROCESSING] 
INFO  13:53:19,550 ProgressMeter -                 |      processed |    time |         per 1M |           |   total | remaining 
INFO  13:53:19,552 ProgressMeter -        Location | active regions | elapsed | active regions | completed | runtime |   runtime 
INFO  13:53:19,552 HaplotypeCaller - Currently, physical phasing is not available when ploidy is different than 2; therefore it won't be performed 
INFO  13:53:19,711 StrandBiasTest - SAM/BAM data was found. Attempting to use read data to calculate strand bias annotations values. 
INFO  13:53:19,711 StrandBiasTest - SAM/BAM data was found. Attempting to use read data to calculate strand bias annotations values. 
Using SSE4.1 accelerated implementation of PairHMM
INFO  13:53:21,527 VectorLoglessPairHMM - libVectorLoglessPairHMM unpacked successfully from GATK jar file 
INFO  13:53:21,528 VectorLoglessPairHMM - Using vectorized implementation of PairHMM 
INFO  13:53:49,556 ProgressMeter - Pf3D7_07_v3:481717              0.0    30.0 s           49.6 w        8.6%     5.8 m       5.3 m 
INFO  13:54:19,558 ProgressMeter - Pf3D7_07_v3:486852              0.0    60.0 s           99.2 w       34.3%     2.9 m     115.0 s 
INFO  13:54:49,559 ProgressMeter - Pf3D7_07_v3:492437              0.0    90.0 s          148.8 w       62.2%     2.4 m      54.0 s 
INFO  13:55:19,561 ProgressMeter - Pf3D7_07_v3:494028              0.0   120.0 s          198.4 w       70.1%     2.9 m      51.0 s 
INFO  13:55:49,563 ProgressMeter - Pf3D7_07_v3:496482              0.0     2.5 m          248.0 w       82.4%     3.0 m      32.0 s 
INFO  13:56:19,564 ProgressMeter - Pf3D7_07_v3:499266              0.0     3.0 m          297.6 w       96.3%     3.1 m       6.0 s 
INFO  13:56:36,889 VectorLoglessPairHMM - Time spent in setup for JNI call : 0.033746580000000005 
INFO  13:56:36,889 PairHMM - Total compute time in PairHMM computeLikelihoods() : 163.053572992 
INFO  13:56:36,890 HaplotypeCaller - Ran local assembly on 0 active regions 
INFO  13:56:36,993 ProgressMeter -            done          20001.0     3.3 m            2.7 h      100.0%     3.3 m       0.0 s 
INFO  13:56:36,994 ProgressMeter - Total runtime 197.44 secs, 3.29 min, 0.05 hours 
INFO  13:56:36,994 MicroScheduler - 7 reads were filtered out during the traversal out of approximately 24416 total reads (0.03%) 
INFO  13:56:36,994 MicroScheduler -   -> 0 reads (0.00% of total) failing BadCigarFilter 
INFO  13:56:36,995 MicroScheduler -   -> 0 reads (0.00% of total) failing DuplicateReadFilter 
INFO  13:56:36,995 MicroScheduler -   -> 0 reads (0.00% of total) failing FailsVendorQualityCheckFilter 
INFO  13:56:36,995 MicroScheduler -   -> 0 reads (0.00% of total) failing HCMappingQualityFilter 
INFO  13:56:36,995 MicroScheduler -   -> 0 reads (0.00% of total) failing MalformedReadFilter 
INFO  13:56:36,996 MicroScheduler -   -> 0 reads (0.00% of total) failing MappingQualityUnavailableFilter 
INFO  13:56:36,997 MicroScheduler -   -> 7 reads (0.03% of total) failing NotPrimaryAlignmentFilter 
INFO  13:56:36,998 MicroScheduler -   -> 0 reads (0.00% of total) failing UnmappedReadFilter 
INFO  13:56:37,663 GATKRunReport - Uploaded run statistics report to AWS S3 
Picked up _JAVA_OPTIONS: -Xmx2000m -Xms50m
INFO  14:31:30,031 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  14:31:30,036 HelpFormatter - The Genome Analysis Toolkit (GATK) v3.4-0-g7e26428, Compiled 2015/05/15 03:25:41 
INFO  14:31:30,036 HelpFormatter - Copyright (c) 2010 The Broad Institute 
INFO  14:31:30,036 HelpFormatter - For support and documentation go to http://www.broadinstitute.org/gatk 
INFO  14:31:30,044 HelpFormatter - Program Args: -T HaplotypeCaller -I bams.3D7.list -ploidy 1 -R /seq/plasmodium/sredmond//refs/Pf3D7_v3.fasta -o haplo.r3D7.ALL.48-50.vcf -gt_mode DISCOVERY -out_mode EMIT_VARIANTS_ONLY -stand_emit_conf 0 -stand_call_conf 30 -L Pf3D7_07_v3:480000-500000 -APO haplo.r3D7.ALL.48-50.apo -globalMAPQ -1 -pcrModel NONE -useFilteredReadsForAnnotations -mmq 0 -nda 
INFO  14:31:30,059 HelpFormatter - Executing as sredmond@a2e-bigmem-0.broadinstitute.org on Linux 2.6.32-504.8.1.el6.x86_64 amd64; Java HotSpot(TM) 64-Bit Server VM 1.8.0-b132. 
INFO  14:31:30,060 HelpFormatter - Date/Time: 2015/07/13 14:31:30 
INFO  14:31:30,061 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  14:31:30,062 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  14:31:30,963 GenomeAnalysisEngine - Strictness is SILENT 
INFO  14:31:31,271 GenomeAnalysisEngine - Downsampling Settings: Method: BY_SAMPLE, Target Coverage: 500 
INFO  14:31:31,284 SAMDataSource$SAMReaders - Initializing SAMRecords in serial 
INFO  14:31:31,358 SAMDataSource$SAMReaders - Done initializing BAM readers: total time 0.07 
INFO  14:31:31,388 IntervalUtils - Processing 20001 bp from intervals 
INFO  14:31:31,523 GenomeAnalysisEngine - Preparing for traversal over 3 BAM files 
INFO  14:31:31,563 GenomeAnalysisEngine - Done preparing for traversal 
INFO  14:31:31,564 ProgressMeter - [INITIALIZATION COMPLETE; STARTING PROCESSING] 
INFO  14:31:31,565 ProgressMeter -                 |      processed |    time |         per 1M |           |   total | remaining 
INFO  14:31:31,565 ProgressMeter -        Location | active regions | elapsed | active regions | completed | runtime |   runtime 
INFO  14:31:31,566 HaplotypeCaller - Currently, physical phasing is not available when ploidy is different than 2; therefore it won't be performed 
INFO  14:31:31,714 StrandBiasTest - SAM/BAM data was found. Attempting to use read data to calculate strand bias annotations values. 
INFO  14:31:31,714 StrandBiasTest - SAM/BAM data was found. Attempting to use read data to calculate strand bias annotations values. 
Using SSE4.1 accelerated implementation of PairHMM
INFO  14:31:34,273 VectorLoglessPairHMM - libVectorLoglessPairHMM unpacked successfully from GATK jar file 
INFO  14:31:34,274 VectorLoglessPairHMM - Using vectorized implementation of PairHMM 
INFO  14:32:01,571 ProgressMeter - Pf3D7_07_v3:481194              0.0    30.0 s           49.6 w        6.0%     8.4 m       7.9 m 
INFO  14:32:31,573 ProgressMeter - Pf3D7_07_v3:484552              0.0    60.0 s           99.2 w       22.8%     4.4 m       3.4 m 
INFO  14:33:01,575 ProgressMeter - Pf3D7_07_v3:488834              0.0    90.0 s          148.8 w       44.2%     3.4 m     113.0 s 
INFO  14:33:31,577 ProgressMeter - Pf3D7_07_v3:492437              0.0   120.0 s          198.4 w       62.2%     3.2 m      72.0 s 
INFO  14:34:01,579 ProgressMeter - Pf3D7_07_v3:493428              0.0     2.5 m          248.0 w       67.1%     3.7 m      73.0 s 
INFO  14:34:31,580 ProgressMeter - Pf3D7_07_v3:495204              0.0     3.0 m          297.6 w       76.0%     3.9 m      56.0 s 
INFO  14:35:01,582 ProgressMeter - Pf3D7_07_v3:496482              0.0     3.5 m          347.3 w       82.4%     4.2 m      44.0 s 
INFO  14:35:31,584 ProgressMeter - Pf3D7_07_v3:498648              0.0     4.0 m          396.9 w       93.2%     4.3 m      17.0 s 
INFO  14:36:01,585 ProgressMeter - Pf3D7_07_v3:499390              0.0     4.5 m          446.5 w       96.9%     4.6 m       8.0 s 
INFO  14:36:22,540 VectorLoglessPairHMM - Time spent in setup for JNI call : 0.040639673 
INFO  14:36:22,541 PairHMM - Total compute time in PairHMM computeLikelihoods() : 243.600002047 
INFO  14:36:22,542 HaplotypeCaller - Ran local assembly on 0 active regions 
INFO  14:36:22,555 ProgressMeter -            done          20001.0     4.8 m            4.0 h      100.0%     4.8 m       0.0 s 
INFO  14:36:22,556 ProgressMeter - Total runtime 290.99 secs, 4.85 min, 0.08 hours 
INFO  14:36:22,557 MicroScheduler - 7 reads were filtered out during the traversal out of approximately 24416 total reads (0.03%) 
INFO  14:36:22,557 MicroScheduler -   -> 0 reads (0.00% of total) failing BadCigarFilter 
INFO  14:36:22,558 MicroScheduler -   -> 0 reads (0.00% of total) failing DuplicateReadFilter 
INFO  14:36:22,558 MicroScheduler -   -> 0 reads (0.00% of total) failing FailsVendorQualityCheckFilter 
INFO  14:36:22,558 MicroScheduler -   -> 0 reads (0.00% of total) failing HCMappingQualityFilter 
INFO  14:36:22,559 MicroScheduler -   -> 0 reads (0.00% of total) failing MalformedReadFilter 
INFO  14:36:22,559 MicroScheduler -   -> 0 reads (0.00% of total) failing MappingQualityUnavailableFilter 
INFO  14:36:22,559 MicroScheduler -   -> 7 reads (0.03% of total) failing NotPrimaryAlignmentFilter 
INFO  14:36:22,560 MicroScheduler -   -> 0 reads (0.00% of total) failing UnmappedReadFilter 
INFO  14:36:23,579 GATKRunReport - Uploaded run statistics report to AWS S3 
Picked up _JAVA_OPTIONS: -Xmx2000m -Xms50m
INFO  16:50:07,946 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  16:50:07,950 HelpFormatter - The Genome Analysis Toolkit (GATK) v3.4-0-g7e26428, Compiled 2015/05/15 03:25:41 
INFO  16:50:07,950 HelpFormatter - Copyright (c) 2010 The Broad Institute 
INFO  16:50:07,951 HelpFormatter - For support and documentation go to http://www.broadinstitute.org/gatk 
INFO  16:50:07,956 HelpFormatter - Program Args: -T HaplotypeCaller -I bams.3D7.list -ploidy 1 -R /seq/plasmodium/sredmond//refs/Pf3D7_v3.fasta -o haplo.r3D7.ALL.vcf -gt_mode DISCOVERY -out_mode EMIT_VARIANTS_ONLY -stand_emit_conf 0 -stand_call_conf 30 -L Pf3D7_07_v3:450000-750000 -APO haplo.r3D7.ALL.apo -globalMAPQ -1 -pcrModel NONE -useFilteredReadsForAnnotations -mmq 0 -nda 
INFO  16:50:07,990 HelpFormatter - Executing as sredmond@a2e-bigmem-0.broadinstitute.org on Linux 2.6.32-504.8.1.el6.x86_64 amd64; Java HotSpot(TM) 64-Bit Server VM 1.8.0-b132. 
INFO  16:50:07,991 HelpFormatter - Date/Time: 2015/07/13 16:50:07 
INFO  16:50:07,992 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  16:50:07,993 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  16:50:08,806 GenomeAnalysisEngine - Strictness is SILENT 
INFO  16:50:08,923 GenomeAnalysisEngine - Downsampling Settings: Method: BY_SAMPLE, Target Coverage: 500 
INFO  16:50:08,935 SAMDataSource$SAMReaders - Initializing SAMRecords in serial 
INFO  16:50:09,608 GATKRunReport - Uploaded run statistics report to AWS S3 
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
##### ERROR MESSAGE: Couldn't read file /seq/plasmodium/sredmond/pfdisco/analyses/fakeNGS/Pf3D7_07_v3:450000-750000_v_Pf3D7.rg.bam because java.io.FileNotFoundException: /seq/plasmodium/sredmond/pfdisco/analyses/fakeNGS/Pf3D7_07_v3:450000-750000_v_Pf3D7.rg.bam (No such file or directory)
##### ERROR ------------------------------------------------------------------------------------------
Picked up _JAVA_OPTIONS: -Xmx2000m -Xms50m
INFO  17:00:59,344 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  17:00:59,347 HelpFormatter - The Genome Analysis Toolkit (GATK) v3.4-0-g7e26428, Compiled 2015/05/15 03:25:41 
INFO  17:00:59,347 HelpFormatter - Copyright (c) 2010 The Broad Institute 
INFO  17:00:59,347 HelpFormatter - For support and documentation go to http://www.broadinstitute.org/gatk 
INFO  17:00:59,352 HelpFormatter - Program Args: -T HaplotypeCaller -I bams.3D7.list -ploidy 1 -R /seq/plasmodium/sredmond//refs/Pf3D7_v3.fasta -o haplo.r3D7.ALL.vcf -gt_mode DISCOVERY -out_mode EMIT_VARIANTS_ONLY -stand_emit_conf 0 -stand_call_conf 30 -L Pf3D7_07_v3:480000-500000 -APO haplo.r3D7.ALL.apo -globalMAPQ -1 -pcrModel NONE -useFilteredReadsForAnnotations -mmq 0 -nda 
INFO  17:00:59,360 HelpFormatter - Executing as sredmond@a2e-bigmem-0.broadinstitute.org on Linux 2.6.32-504.8.1.el6.x86_64 amd64; Java HotSpot(TM) 64-Bit Server VM 1.8.0-b132. 
INFO  17:00:59,361 HelpFormatter - Date/Time: 2015/07/13 17:00:59 
INFO  17:00:59,361 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  17:00:59,363 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  17:01:00,026 GenomeAnalysisEngine - Strictness is SILENT 
INFO  17:01:00,261 GenomeAnalysisEngine - Downsampling Settings: Method: BY_SAMPLE, Target Coverage: 500 
INFO  17:01:00,271 SAMDataSource$SAMReaders - Initializing SAMRecords in serial 
INFO  17:01:00,342 SAMDataSource$SAMReaders - Done initializing BAM readers: total time 0.07 
INFO  17:01:00,364 IntervalUtils - Processing 20001 bp from intervals 
INFO  17:01:00,467 GenomeAnalysisEngine - Preparing for traversal over 3 BAM files 
INFO  17:01:00,496 GenomeAnalysisEngine - Done preparing for traversal 
INFO  17:01:00,496 ProgressMeter - [INITIALIZATION COMPLETE; STARTING PROCESSING] 
INFO  17:01:00,497 ProgressMeter -                 |      processed |    time |         per 1M |           |   total | remaining 
INFO  17:01:00,497 ProgressMeter -        Location | active regions | elapsed | active regions | completed | runtime |   runtime 
INFO  17:01:00,497 HaplotypeCaller - Currently, physical phasing is not available when ploidy is different than 2; therefore it won't be performed 
INFO  17:01:00,596 StrandBiasTest - SAM/BAM data was found. Attempting to use read data to calculate strand bias annotations values. 
INFO  17:01:00,597 StrandBiasTest - SAM/BAM data was found. Attempting to use read data to calculate strand bias annotations values. 
Using SSE4.1 accelerated implementation of PairHMM
INFO  17:01:02,175 VectorLoglessPairHMM - libVectorLoglessPairHMM unpacked successfully from GATK jar file 
INFO  17:01:02,176 VectorLoglessPairHMM - Using vectorized implementation of PairHMM 
WARN  17:01:05,481 HaplotypeScore - Annotation will not be calculated, must be called from UnifiedGenotyper 
WARN  17:01:05,482 InbreedingCoeff - Annotation will not be calculated, must provide a valid PED file (-ped) from the command line. 
INFO  17:01:30,502 ProgressMeter - Pf3D7_07_v3:481905              0.0    30.0 s           49.6 w        9.5%     5.2 m       4.7 m 
WARN  17:01:42,102 AnnotationUtils - Annotation will not be calculated, genotype is not called 
INFO  17:02:00,504 ProgressMeter - Pf3D7_07_v3:486852              0.0    60.0 s           99.2 w       34.3%     2.9 m     115.0 s 
INFO  17:02:30,505 ProgressMeter - Pf3D7_07_v3:492603              0.0    90.0 s          148.8 w       63.0%     2.4 m      52.0 s 
INFO  17:03:00,507 ProgressMeter - Pf3D7_07_v3:494304              0.0   120.0 s          198.4 w       71.5%     2.8 m      47.0 s 
Picked up _JAVA_OPTIONS: -Xmx2000m -Xms50m
INFO  17:05:19,440 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  17:05:19,443 HelpFormatter - The Genome Analysis Toolkit (GATK) v3.4-0-g7e26428, Compiled 2015/05/15 03:25:41 
INFO  17:05:19,443 HelpFormatter - Copyright (c) 2010 The Broad Institute 
INFO  17:05:19,443 HelpFormatter - For support and documentation go to http://www.broadinstitute.org/gatk 
INFO  17:05:19,448 HelpFormatter - Program Args: -T HaplotypeCaller -I bams.3D7.list -ploidy 1 -R /seq/plasmodium/sredmond//refs/Pf3D7_v3.fasta -o haplo.r3D7.ALL.vcf -gt_mode DISCOVERY -out_mode EMIT_VARIANTS_ONLY -stand_emit_conf 0 -stand_call_conf 30 -L Pf3D7_07_v3:480000-500000 -APO haplo.r3D7.ALL.apo -globalMAPQ -1 -pcrModel NONE -useFilteredReadsForAnnotations -mmq 0 -nda 
INFO  17:05:19,464 HelpFormatter - Executing as sredmond@a2e-bigmem-0.broadinstitute.org on Linux 2.6.32-504.8.1.el6.x86_64 amd64; Java HotSpot(TM) 64-Bit Server VM 1.8.0-b132. 
INFO  17:05:19,465 HelpFormatter - Date/Time: 2015/07/13 17:05:19 
INFO  17:05:19,465 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  17:05:19,466 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  17:05:20,151 GenomeAnalysisEngine - Strictness is SILENT 
INFO  17:05:20,461 GenomeAnalysisEngine - Downsampling Settings: Method: BY_SAMPLE, Target Coverage: 500 
INFO  17:05:20,470 SAMDataSource$SAMReaders - Initializing SAMRecords in serial 
INFO  17:05:20,559 SAMDataSource$SAMReaders - Done initializing BAM readers: total time 0.09 
INFO  17:05:20,588 IntervalUtils - Processing 20001 bp from intervals 
INFO  17:05:20,703 GenomeAnalysisEngine - Preparing for traversal over 3 BAM files 
INFO  17:05:20,812 GenomeAnalysisEngine - Done preparing for traversal 
INFO  17:05:20,813 ProgressMeter - [INITIALIZATION COMPLETE; STARTING PROCESSING] 
INFO  17:05:20,814 ProgressMeter -                 |      processed |    time |         per 1M |           |   total | remaining 
INFO  17:05:20,815 ProgressMeter -        Location | active regions | elapsed | active regions | completed | runtime |   runtime 
INFO  17:05:20,815 HaplotypeCaller - Currently, physical phasing is not available when ploidy is different than 2; therefore it won't be performed 
INFO  17:05:21,008 StrandBiasTest - SAM/BAM data was found. Attempting to use read data to calculate strand bias annotations values. 
INFO  17:05:21,008 StrandBiasTest - SAM/BAM data was found. Attempting to use read data to calculate strand bias annotations values. 
Using SSE4.1 accelerated implementation of PairHMM
INFO  17:05:22,733 VectorLoglessPairHMM - libVectorLoglessPairHMM unpacked successfully from GATK jar file 
INFO  17:05:22,734 VectorLoglessPairHMM - Using vectorized implementation of PairHMM 
WARN  17:05:26,071 HaplotypeScore - Annotation will not be calculated, must be called from UnifiedGenotyper 
WARN  17:05:26,072 InbreedingCoeff - Annotation will not be calculated, must provide a valid PED file (-ped) from the command line. 
INFO  17:05:50,819 ProgressMeter - Pf3D7_07_v3:481905              0.0    30.0 s           49.6 w        9.5%     5.2 m       4.7 m 
WARN  17:06:03,356 AnnotationUtils - Annotation will not be calculated, genotype is not called 
INFO  17:06:20,832 ProgressMeter - Pf3D7_07_v3:486852              0.0    60.0 s           99.2 w       34.3%     2.9 m     115.0 s 
INFO  17:06:50,834 ProgressMeter - Pf3D7_07_v3:492437              0.0    90.0 s          148.8 w       62.2%     2.4 m      54.0 s 
INFO  17:07:20,836 ProgressMeter - Pf3D7_07_v3:494021              0.0   120.0 s          198.5 w       70.1%     2.9 m      51.0 s 
INFO  17:07:50,838 ProgressMeter - Pf3D7_07_v3:496482              0.0     2.5 m          248.1 w       82.4%     3.0 m      32.0 s 
INFO  17:08:20,840 ProgressMeter - Pf3D7_07_v3:498648              0.0     3.0 m          297.7 w       93.2%     3.2 m      13.0 s 
INFO  17:08:45,473 VectorLoglessPairHMM - Time spent in setup for JNI call : 0.033151294000000005 
INFO  17:08:45,474 PairHMM - Total compute time in PairHMM computeLikelihoods() : 164.85849562 
INFO  17:08:45,474 HaplotypeCaller - Ran local assembly on 58 active regions 
INFO  17:08:45,488 ProgressMeter -            done          20001.0     3.4 m            2.8 h      100.0%     3.4 m       0.0 s 
INFO  17:08:45,488 ProgressMeter - Total runtime 204.68 secs, 3.41 min, 0.06 hours 
INFO  17:08:45,489 MicroScheduler - 8 reads were filtered out during the traversal out of approximately 24472 total reads (0.03%) 
INFO  17:08:45,489 MicroScheduler -   -> 0 reads (0.00% of total) failing BadCigarFilter 
INFO  17:08:45,489 MicroScheduler -   -> 0 reads (0.00% of total) failing DuplicateReadFilter 
INFO  17:08:45,490 MicroScheduler -   -> 0 reads (0.00% of total) failing FailsVendorQualityCheckFilter 
INFO  17:08:45,490 MicroScheduler -   -> 0 reads (0.00% of total) failing HCMappingQualityFilter 
INFO  17:08:45,490 MicroScheduler -   -> 0 reads (0.00% of total) failing MalformedReadFilter 
INFO  17:08:45,491 MicroScheduler -   -> 0 reads (0.00% of total) failing MappingQualityUnavailableFilter 
INFO  17:08:45,491 MicroScheduler -   -> 8 reads (0.03% of total) failing NotPrimaryAlignmentFilter 
INFO  17:08:45,491 MicroScheduler -   -> 0 reads (0.00% of total) failing UnmappedReadFilter 
INFO  17:08:46,363 GATKRunReport - Uploaded run statistics report to AWS S3 
Picked up _JAVA_OPTIONS: -Xmx2000m -Xms50m
INFO  17:26:33,683 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  17:26:33,686 HelpFormatter - The Genome Analysis Toolkit (GATK) v3.4-0-g7e26428, Compiled 2015/05/15 03:25:41 
INFO  17:26:33,686 HelpFormatter - Copyright (c) 2010 The Broad Institute 
INFO  17:26:33,686 HelpFormatter - For support and documentation go to http://www.broadinstitute.org/gatk 
INFO  17:26:33,691 HelpFormatter - Program Args: -T HaplotypeCaller -I bams.3D7.list -ploidy 1 -R /seq/plasmodium/sredmond//refs/Pf3D7_v3.fasta -o haplo.r3D7.ALL.vcf -gt_mode DISCOVERY -out_mode EMIT_VARIANTS_ONLY -stand_emit_conf 0 -stand_call_conf 30 -L Pf3D7_07_v3:450000-750000 -APO haplo.r3D7.ALL.apo -globalMAPQ -1 -pcrModel NONE -useFilteredReadsForAnnotations -mmq 0 -nda 
INFO  17:26:33,700 HelpFormatter - Executing as sredmond@a2e-bigmem-0.broadinstitute.org on Linux 2.6.32-504.8.1.el6.x86_64 amd64; Java HotSpot(TM) 64-Bit Server VM 1.8.0-b132. 
INFO  17:26:33,701 HelpFormatter - Date/Time: 2015/07/13 17:26:33 
INFO  17:26:33,701 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  17:26:33,701 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  17:26:34,613 GenomeAnalysisEngine - Strictness is SILENT 
INFO  17:26:34,765 GenomeAnalysisEngine - Downsampling Settings: Method: BY_SAMPLE, Target Coverage: 500 
INFO  17:26:34,775 SAMDataSource$SAMReaders - Initializing SAMRecords in serial 
INFO  17:26:34,867 SAMDataSource$SAMReaders - Done initializing BAM readers: total time 0.09 
INFO  17:26:34,893 IntervalUtils - Processing 300001 bp from intervals 
INFO  17:26:35,014 GenomeAnalysisEngine - Preparing for traversal over 3 BAM files 
INFO  17:26:35,056 GenomeAnalysisEngine - Done preparing for traversal 
INFO  17:26:35,057 ProgressMeter - [INITIALIZATION COMPLETE; STARTING PROCESSING] 
INFO  17:26:35,058 ProgressMeter -                 |      processed |    time |         per 1M |           |   total | remaining 
INFO  17:26:35,058 ProgressMeter -        Location | active regions | elapsed | active regions | completed | runtime |   runtime 
INFO  17:26:35,059 HaplotypeCaller - Currently, physical phasing is not available when ploidy is different than 2; therefore it won't be performed 
INFO  17:26:35,209 StrandBiasTest - SAM/BAM data was found. Attempting to use read data to calculate strand bias annotations values. 
INFO  17:26:35,210 StrandBiasTest - SAM/BAM data was found. Attempting to use read data to calculate strand bias annotations values. 
Using SSE4.1 accelerated implementation of PairHMM
INFO  17:26:36,799 VectorLoglessPairHMM - libVectorLoglessPairHMM unpacked successfully from GATK jar file 
INFO  17:26:36,800 VectorLoglessPairHMM - Using vectorized implementation of PairHMM 
WARN  17:26:37,760 HaplotypeScore - Annotation will not be calculated, must be called from UnifiedGenotyper 
WARN  17:26:37,760 InbreedingCoeff - Annotation will not be calculated, must provide a valid PED file (-ped) from the command line. 
WARN  17:26:37,763 AnnotationUtils - Annotation will not be calculated, genotype is not called 
INFO  17:27:05,063 ProgressMeter - Pf3D7_07_v3:454987              0.0    30.0 s           49.6 w        1.7%    30.1 m      29.6 m 
INFO  17:27:35,065 ProgressMeter - Pf3D7_07_v3:459474              0.0    60.0 s           99.2 w        3.2%    31.7 m      30.7 m 
INFO  17:28:05,067 ProgressMeter - Pf3D7_07_v3:463227              0.0    90.0 s          148.8 w        4.4%    34.0 m      32.5 m 
INFO  17:28:35,068 ProgressMeter - Pf3D7_07_v3:465649              0.0   120.0 s          198.4 w        5.2%    38.3 m      36.3 m 
INFO  17:29:05,070 ProgressMeter - Pf3D7_07_v3:467011              0.0     2.5 m          248.0 w        5.7%    44.1 m      41.6 m 
INFO  17:29:35,071 ProgressMeter - Pf3D7_07_v3:469238              0.0     3.0 m          297.6 w        6.4%    46.8 m      43.8 m 
INFO  17:30:05,072 ProgressMeter - Pf3D7_07_v3:473268              0.0     3.5 m          347.2 w        7.8%    45.1 m      41.6 m 
INFO  17:30:35,074 ProgressMeter - Pf3D7_07_v3:474755              0.0     4.0 m          396.9 w        8.3%    48.5 m      44.5 m 
INFO  17:31:05,076 ProgressMeter - Pf3D7_07_v3:476554              0.0     4.5 m          446.5 w        8.9%    50.8 m      46.3 m 
INFO  17:31:35,078 ProgressMeter - Pf3D7_07_v3:478120              0.0     5.0 m          496.1 w        9.4%    53.3 m      48.3 m 
INFO  17:32:05,079 ProgressMeter - Pf3D7_07_v3:480508              0.0     5.5 m          545.7 w       10.2%    54.1 m      48.6 m 
INFO  17:32:35,080 ProgressMeter - Pf3D7_07_v3:481905              0.0     6.0 m          595.3 w       10.6%    56.4 m      50.4 m 
INFO  17:33:05,082 ProgressMeter - Pf3D7_07_v3:484745              0.0     6.5 m          644.9 w       11.6%    56.1 m      49.6 m 
INFO  17:33:35,083 ProgressMeter - Pf3D7_07_v3:490408              0.0     7.0 m          694.5 w       13.5%    52.0 m      45.0 m 
INFO  17:34:05,084 ProgressMeter - Pf3D7_07_v3:492603              0.0     7.5 m          744.1 w       14.2%    52.8 m      45.3 m 
INFO  17:34:35,086 ProgressMeter - Pf3D7_07_v3:493856              0.0     8.0 m          793.7 w       14.6%    54.7 m      46.7 m 
INFO  17:35:05,087 ProgressMeter - Pf3D7_07_v3:495197              0.0     8.5 m          843.3 w       15.1%    56.4 m      47.9 m 
INFO  17:35:35,088 ProgressMeter - Pf3D7_07_v3:496482              0.0     9.0 m          892.9 w       15.5%    58.1 m      49.1 m 
INFO  17:36:05,089 ProgressMeter - Pf3D7_07_v3:498648              0.0     9.5 m          942.5 w       16.2%    58.6 m      49.1 m 
INFO  17:36:35,091 ProgressMeter - Pf3D7_07_v3:502949              0.0    10.0 m          992.1 w       17.6%    56.7 m      46.7 m 
WARN  17:36:46,447 ExactAFCalculator - this tool is currently set to genotype at most 6 alternate alleles in a given context, but the context at Pf3D7_07_v3:503153 has 9 alternate alleles so only the top alleles will be used; see the --max_alternate_alleles argument 
INFO  17:37:05,092 ProgressMeter - Pf3D7_07_v3:504462              0.0    10.5 m         1041.7 w       18.2%    57.8 m      47.3 m 
INFO  17:37:35,094 ProgressMeter - Pf3D7_07_v3:505401              0.0    11.0 m         1091.3 w       18.5%    59.6 m      48.6 m 
INFO  17:38:05,095 ProgressMeter - Pf3D7_07_v3:508269              0.0    11.5 m         1140.9 w       19.4%    59.2 m      47.7 m 
INFO  17:38:35,096 ProgressMeter - Pf3D7_07_v3:509948              0.0    12.0 m         1190.5 w       20.0%    60.1 m      48.1 m 
INFO  17:39:05,098 ProgressMeter - Pf3D7_07_v3:511210              0.0    12.5 m         1240.1 w       20.4%    61.3 m      48.8 m 
INFO  17:39:35,099 ProgressMeter - Pf3D7_07_v3:515680              0.0    13.0 m         1289.8 w       21.9%    59.4 m      46.4 m 
INFO  17:40:05,100 ProgressMeter - Pf3D7_07_v3:521247              0.0    13.5 m         1339.4 w       23.7%    56.8 m      43.3 m 
INFO  17:40:35,101 ProgressMeter - Pf3D7_07_v3:529416              0.0    14.0 m         1389.0 w       26.5%    52.9 m      38.9 m 
INFO  17:41:05,102 ProgressMeter - Pf3D7_07_v3:533793              0.0    14.5 m         1438.6 w       27.9%    51.9 m      37.4 m 
INFO  17:41:35,104 ProgressMeter - Pf3D7_07_v3:547567              0.0    15.0 m         1488.2 w       32.5%    46.1 m      31.1 m 
INFO  17:42:05,105 ProgressMeter - Pf3D7_07_v3:551488              0.0    15.5 m         1537.8 w       33.8%    45.8 m      30.3 m 
INFO  17:42:35,107 ProgressMeter - Pf3D7_07_v3:552645              0.0    16.0 m         1587.4 w       34.2%    46.8 m      30.8 m 
INFO  17:43:05,108 ProgressMeter - Pf3D7_07_v3:554360              0.0    16.5 m         1637.0 w       34.8%    47.4 m      30.9 m 
INFO  17:43:35,109 ProgressMeter - Pf3D7_07_v3:560618              0.0    17.0 m         1686.6 w       36.9%    46.1 m      29.1 m 
INFO  17:44:05,111 ProgressMeter - Pf3D7_07_v3:564395              0.0    17.5 m         1736.2 w       38.1%    45.9 m      28.4 m 
INFO  17:44:35,112 ProgressMeter - Pf3D7_07_v3:566664              0.0    18.0 m         1785.8 w       38.9%    46.3 m      28.3 m 
INFO  17:45:05,113 ProgressMeter - Pf3D7_07_v3:569166              0.0    18.5 m         1835.4 w       39.7%    46.6 m      28.1 m 
INFO  17:45:35,114 ProgressMeter - Pf3D7_07_v3:572795              0.0    19.0 m         1885.0 w       40.9%    46.4 m      27.4 m 
INFO  17:46:05,116 ProgressMeter - Pf3D7_07_v3:574097              0.0    19.5 m         1934.6 w       41.4%    47.1 m      27.6 m 
INFO  17:46:35,117 ProgressMeter - Pf3D7_07_v3:576493              0.0    20.0 m         1984.2 w       42.2%    47.4 m      27.4 m 
INFO  17:47:05,118 ProgressMeter - Pf3D7_07_v3:579756              0.0    20.5 m         2033.8 w       43.3%    47.4 m      26.9 m 
INFO  17:47:35,119 ProgressMeter - Pf3D7_07_v3:581542              0.0    21.0 m         2083.4 w       43.8%    47.9 m      26.9 m 
INFO  17:48:05,121 ProgressMeter - Pf3D7_07_v3:585588              0.0    21.5 m         2133.0 w       45.2%    47.6 m      26.1 m 
INFO  17:48:35,122 ProgressMeter - Pf3D7_07_v3:590544              0.0    22.0 m         2182.6 w       46.8%    47.0 m      25.0 m 
INFO  17:49:05,123 ProgressMeter - Pf3D7_07_v3:595492              0.0    22.5 m         2232.3 w       48.5%    46.4 m      23.9 m 
INFO  17:49:35,124 ProgressMeter - Pf3D7_07_v3:599139              0.0    23.0 m         2281.9 w       49.7%    46.3 m      23.3 m 
INFO  17:50:05,125 ProgressMeter - Pf3D7_07_v3:599851              0.0    23.5 m         2331.5 w       50.0%    47.0 m      23.5 m 
INFO  17:50:35,126 ProgressMeter - Pf3D7_07_v3:600532              0.0    24.0 m         2381.1 w       50.2%    47.8 m      23.8 m 
INFO  17:51:05,128 ProgressMeter - Pf3D7_07_v3:600986              0.0    24.5 m         2430.7 w       50.3%    48.7 m      24.2 m 
INFO  17:51:35,129 ProgressMeter - Pf3D7_07_v3:601894              0.0    25.0 m         2480.3 w       50.6%    49.4 m      24.4 m 
INFO  17:52:05,130 ProgressMeter - Pf3D7_07_v3:602687              0.0    25.5 m         2529.9 w       50.9%    50.1 m      24.6 m 
INFO  17:52:35,132 ProgressMeter - Pf3D7_07_v3:604167              0.0    26.0 m         2579.5 w       51.4%    50.6 m      24.6 m 
INFO  17:53:05,133 ProgressMeter - Pf3D7_07_v3:604516              0.0    26.5 m         2629.1 w       51.5%    51.5 m      25.0 m 
INFO  17:53:35,134 ProgressMeter - Pf3D7_07_v3:605175              0.0    27.0 m         2678.7 w       51.7%    52.2 m      25.2 m 
INFO  17:54:05,135 ProgressMeter - Pf3D7_07_v3:605956              0.0    27.5 m         2728.3 w       52.0%    52.9 m      25.4 m 
INFO  17:54:35,136 ProgressMeter - Pf3D7_07_v3:606912              0.0    28.0 m         2777.9 w       52.3%    53.5 m      25.5 m 
INFO  17:55:05,137 ProgressMeter - Pf3D7_07_v3:608993              0.0    28.5 m         2827.5 w       53.0%    53.8 m      25.3 m 
INFO  17:55:35,141 ProgressMeter - Pf3D7_07_v3:610329              0.0    29.0 m         2877.1 w       53.4%    54.3 m      25.3 m 
INFO  17:56:05,142 ProgressMeter - Pf3D7_07_v3:611896              0.0    29.5 m         2926.7 w       54.0%    54.7 m      25.2 m 
INFO  17:56:35,143 ProgressMeter - Pf3D7_07_v3:613269              0.0    30.0 m         2976.3 w       54.4%    55.1 m      25.1 m 
INFO  17:57:05,145 ProgressMeter - Pf3D7_07_v3:614880              0.0    30.5 m         3025.9 w       55.0%    55.5 m      25.0 m 
INFO  17:57:35,146 ProgressMeter - Pf3D7_07_v3:617662              0.0    31.0 m         3075.5 w       55.9%    55.5 m      24.5 m 
INFO  17:58:05,147 ProgressMeter - Pf3D7_07_v3:619642              0.0    31.5 m         3125.1 w       56.5%    55.7 m      24.2 m 
INFO  17:58:35,163 ProgressMeter - Pf3D7_07_v3:621429              0.0    32.0 m         3174.8 w       57.1%    56.0 m      24.0 m 
INFO  17:59:05,165 ProgressMeter - Pf3D7_07_v3:625963              0.0    32.5 m         3224.4 w       58.7%    55.4 m      22.9 m 
INFO  17:59:35,166 ProgressMeter - Pf3D7_07_v3:627505              0.0    33.0 m         3274.0 w       59.2%    55.8 m      22.8 m 
INFO  18:00:05,167 ProgressMeter - Pf3D7_07_v3:628626              0.0    33.5 m         3323.6 w       59.5%    56.3 m      22.8 m 
INFO  18:00:35,168 ProgressMeter - Pf3D7_07_v3:634679              0.0    34.0 m         3373.2 w       61.6%    55.2 m      21.2 m 
INFO  18:01:05,169 ProgressMeter - Pf3D7_07_v3:636334              0.0    34.5 m         3422.8 w       62.1%    55.5 m      21.0 m 
WARN  18:01:16,440 ExactAFCalculator - this tool is currently set to genotype at most 6 alternate alleles in a given context, but the context at Pf3D7_07_v3:636498 has 8 alternate alleles so only the top alleles will be used; see the --max_alternate_alleles argument 
INFO  18:01:35,181 ProgressMeter - Pf3D7_07_v3:639212              0.0    35.0 m         3472.4 w       63.1%    55.5 m      20.5 m 
INFO  18:02:05,182 ProgressMeter - Pf3D7_07_v3:640496              0.0    35.5 m         3522.0 w       63.5%    55.9 m      20.4 m 
INFO  18:02:35,183 ProgressMeter - Pf3D7_07_v3:644343              0.0    36.0 m         3571.6 w       64.8%    55.6 m      19.6 m 
INFO  18:03:05,184 ProgressMeter - Pf3D7_07_v3:646166              0.0    36.5 m         3621.2 w       65.4%    55.8 m      19.3 m 
INFO  18:03:35,185 ProgressMeter - Pf3D7_07_v3:649909              0.0    37.0 m         3670.8 w       66.6%    55.5 m      18.5 m 
INFO  18:04:05,186 ProgressMeter - Pf3D7_07_v3:652486              0.0    37.5 m         3720.5 w       67.5%    55.6 m      18.1 m 
INFO  18:04:35,187 ProgressMeter - Pf3D7_07_v3:656774              0.0    38.0 m         3770.1 w       68.9%    55.1 m      17.1 m 
INFO  18:05:05,188 ProgressMeter - Pf3D7_07_v3:664879              0.0    38.5 m         3819.7 w       71.6%    53.8 m      15.3 m 
INFO  18:05:35,189 ProgressMeter - Pf3D7_07_v3:666771              0.0    39.0 m         3869.3 w       72.3%    54.0 m      15.0 m 
INFO  18:06:05,191 ProgressMeter - Pf3D7_07_v3:668212              0.0    39.5 m         3918.9 w       72.7%    54.3 m      14.8 m 
INFO  18:06:35,192 ProgressMeter - Pf3D7_07_v3:670670              0.0    40.0 m         3968.5 w       73.6%    54.4 m      14.4 m 
INFO  18:07:05,193 ProgressMeter - Pf3D7_07_v3:674701              0.0    40.5 m         4018.1 w       74.9%    54.1 m      13.6 m 
INFO  18:07:35,195 ProgressMeter - Pf3D7_07_v3:675610              0.0    41.0 m         4067.7 w       75.2%    54.5 m      13.5 m 
INFO  18:08:05,196 ProgressMeter - Pf3D7_07_v3:677276              0.0    41.5 m         4117.3 w       75.8%    54.8 m      13.3 m 
INFO  18:08:35,197 ProgressMeter - Pf3D7_07_v3:680837              0.0    42.0 m         4166.9 w       76.9%    54.6 m      12.6 m 
INFO  18:09:05,198 ProgressMeter - Pf3D7_07_v3:683613              0.0    42.5 m         4216.5 w       77.9%    54.6 m      12.1 m 
INFO  18:09:35,223 ProgressMeter - Pf3D7_07_v3:685039              0.0    43.0 m         4266.1 w       78.3%    54.9 m      11.9 m 
WARN  18:09:52,747 ExactAFCalculator - this tool is currently set to genotype at most 6 alternate alleles in a given context, but the context at Pf3D7_07_v3:685341 has 7 alternate alleles so only the top alleles will be used; see the --max_alternate_alleles argument 
INFO  18:10:05,224 ProgressMeter - Pf3D7_07_v3:687352              0.0    43.5 m         4315.8 w       79.1%    55.0 m      11.5 m 
INFO  18:10:35,225 ProgressMeter - Pf3D7_07_v3:689390              0.0    44.0 m         4365.4 w       79.8%    55.1 m      11.1 m 
INFO  18:11:05,226 ProgressMeter - Pf3D7_07_v3:693669              0.0    44.5 m         4415.0 w       81.2%    54.8 m      10.3 m 
INFO  18:11:35,227 ProgressMeter - Pf3D7_07_v3:694643              0.0    45.0 m         4464.6 w       81.5%    55.2 m      10.2 m 
INFO  18:12:05,229 ProgressMeter - Pf3D7_07_v3:695905              0.0    45.5 m         4514.2 w       82.0%    55.5 m      10.0 m 
INFO  18:12:35,230 ProgressMeter - Pf3D7_07_v3:697971              0.0    46.0 m         4563.8 w       82.7%    55.7 m       9.7 m 
INFO  18:13:05,231 ProgressMeter - Pf3D7_07_v3:701264              0.0    46.5 m         4613.4 w       83.8%    55.5 m       9.0 m 
INFO  18:13:35,232 ProgressMeter - Pf3D7_07_v3:703295              0.0    47.0 m         4663.0 w       84.4%    55.7 m       8.7 m 
INFO  18:14:05,233 ProgressMeter - Pf3D7_07_v3:706436              0.0    47.5 m         4712.6 w       85.5%    55.6 m       8.1 m 
INFO  18:14:35,235 ProgressMeter - Pf3D7_07_v3:706492              0.0    48.0 m         4762.2 w       85.5%    56.1 m       8.1 m 
INFO  18:15:05,236 ProgressMeter - Pf3D7_07_v3:707978              0.0    48.5 m         4811.8 w       86.0%    56.4 m       7.9 m 
INFO  18:15:35,238 ProgressMeter - Pf3D7_07_v3:711815              0.0    49.0 m         4861.4 w       87.3%    56.1 m       7.1 m 
INFO  18:16:05,240 ProgressMeter - Pf3D7_07_v3:713478              0.0    49.5 m         4911.0 w       87.8%    56.4 m       6.9 m 
INFO  18:16:35,242 ProgressMeter - Pf3D7_07_v3:716311              0.0    50.0 m         4960.6 w       88.8%    56.3 m       6.3 m 
INFO  18:17:05,243 ProgressMeter - Pf3D7_07_v3:717112              0.0    50.5 m         5010.2 w       89.0%    56.7 m       6.2 m 
INFO  18:17:35,244 ProgressMeter - Pf3D7_07_v3:719118              0.0    51.0 m         5059.8 w       89.7%    56.9 m       5.9 m 
INFO  18:18:05,245 ProgressMeter - Pf3D7_07_v3:719797              0.0    51.5 m         5109.4 w       89.9%    57.3 m       5.8 m 
INFO  18:18:35,246 ProgressMeter - Pf3D7_07_v3:721926              0.0    52.0 m         5159.0 w       90.6%    57.4 m       5.4 m 
INFO  18:19:05,247 ProgressMeter - Pf3D7_07_v3:723696              0.0    52.5 m         5208.6 w       91.2%    57.5 m       5.0 m 
INFO  18:19:35,248 ProgressMeter - Pf3D7_07_v3:724779              0.0    53.0 m         5258.3 w       91.6%    57.9 m       4.9 m 
INFO  18:20:05,249 ProgressMeter - Pf3D7_07_v3:726831              0.0    53.5 m         5307.9 w       92.3%    58.0 m       4.5 m 
INFO  18:20:35,250 ProgressMeter - Pf3D7_07_v3:727948              0.0    54.0 m         5357.5 w       92.6%    58.3 m       4.3 m 
INFO  18:21:05,251 ProgressMeter - Pf3D7_07_v3:729913              0.0    54.5 m         5407.1 w       93.3%    58.4 m       3.9 m 
INFO  18:21:35,253 ProgressMeter - Pf3D7_07_v3:732260              0.0    55.0 m         5456.7 w       94.1%    58.5 m       3.5 m 
INFO  18:22:05,254 ProgressMeter - Pf3D7_07_v3:733951              0.0    55.5 m         5506.3 w       94.7%    58.6 m       3.1 m 
INFO  18:22:35,255 ProgressMeter - Pf3D7_07_v3:735237              0.0    56.0 m         5555.9 w       95.1%    58.9 m       2.9 m 
INFO  18:23:05,256 ProgressMeter - Pf3D7_07_v3:736064              0.0    56.5 m         5605.5 w       95.4%    59.3 m       2.8 m 
INFO  18:23:35,257 ProgressMeter - Pf3D7_07_v3:738269              0.0    57.0 m         5655.1 w       96.1%    59.3 m       2.3 m 
INFO  18:24:05,258 ProgressMeter - Pf3D7_07_v3:739912              0.0    57.5 m         5704.7 w       96.6%    59.5 m     120.0 s 
WARN  18:24:24,490 ExactAFCalculator - this tool is currently set to genotype at most 6 alternate alleles in a given context, but the context at Pf3D7_07_v3:740363 has 8 alternate alleles so only the top alleles will be used; see the --max_alternate_alleles argument 
INFO  18:24:35,259 ProgressMeter - Pf3D7_07_v3:741075              0.0    58.0 m         5754.3 w       97.0%    59.8 m     106.0 s 
INFO  18:25:05,260 ProgressMeter - Pf3D7_07_v3:741867              0.0    58.5 m         5803.9 w       97.3%    60.1 m      97.0 s 
INFO  18:25:35,263 ProgressMeter - Pf3D7_07_v3:743192              0.0    59.0 m         5853.5 w       97.7%    60.4 m      82.0 s 
INFO  18:26:05,280 ProgressMeter - Pf3D7_07_v3:745754              0.0    59.5 m         5903.1 w       98.6%    60.4 m      51.0 s 
INFO  18:26:35,282 ProgressMeter - Pf3D7_07_v3:747579              0.0    60.0 m         5952.8 w       99.2%    60.5 m      29.0 s 
INFO  18:27:05,283 ProgressMeter - Pf3D7_07_v3:748432              0.0    60.5 m         6002.4 w       99.5%    60.8 m      19.0 s 
INFO  18:27:23,144 VectorLoglessPairHMM - Time spent in setup for JNI call : 0.425720686 
INFO  18:27:23,145 PairHMM - Total compute time in PairHMM computeLikelihoods() : 2984.396098953 
INFO  18:27:23,146 HaplotypeCaller - Ran local assembly on 635 active regions 
INFO  18:27:23,172 ProgressMeter -            done         300001.0    60.8 m            3.4 h      100.0%    60.8 m       0.0 s 
INFO  18:27:23,173 ProgressMeter - Total runtime 3648.12 secs, 60.80 min, 1.01 hours 
INFO  18:27:23,174 MicroScheduler - 2168 reads were filtered out during the traversal out of approximately 289986 total reads (0.75%) 
INFO  18:27:23,174 MicroScheduler -   -> 0 reads (0.00% of total) failing BadCigarFilter 
INFO  18:27:23,175 MicroScheduler -   -> 0 reads (0.00% of total) failing DuplicateReadFilter 
INFO  18:27:23,175 MicroScheduler -   -> 0 reads (0.00% of total) failing FailsVendorQualityCheckFilter 
INFO  18:27:23,175 MicroScheduler -   -> 0 reads (0.00% of total) failing HCMappingQualityFilter 
INFO  18:27:23,175 MicroScheduler -   -> 0 reads (0.00% of total) failing MalformedReadFilter 
INFO  18:27:23,176 MicroScheduler -   -> 0 reads (0.00% of total) failing MappingQualityUnavailableFilter 
INFO  18:27:23,176 MicroScheduler -   -> 2019 reads (0.70% of total) failing NotPrimaryAlignmentFilter 
INFO  18:27:23,176 MicroScheduler -   -> 149 reads (0.05% of total) failing UnmappedReadFilter 
INFO  18:27:24,105 GATKRunReport - Uploaded run statistics report to AWS S3 
Picked up _JAVA_OPTIONS: -Xmx2000m -Xms50m
INFO  19:44:38,307 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  19:44:38,311 HelpFormatter - The Genome Analysis Toolkit (GATK) v3.4-0-g7e26428, Compiled 2015/05/15 03:25:41 
INFO  19:44:38,312 HelpFormatter - Copyright (c) 2010 The Broad Institute 
INFO  19:44:38,312 HelpFormatter - For support and documentation go to http://www.broadinstitute.org/gatk 
INFO  19:44:38,319 HelpFormatter - Program Args: -T HaplotypeCaller -I bams.DD2.list -ploidy 1 -R /seq/plasmodium/sredmond//refs/PfDD2_v1.fasta -o haplo.rDD2.ALL.vcf -gt_mode DISCOVERY -out_mode EMIT_VARIANTS_ONLY -stand_emit_conf 0 -stand_call_conf 30 -L PfDD2_07_v1:450000-750000 -APO haplo.rDD2.ALL.apo -globalMAPQ -1 -pcrModel NONE -useFilteredReadsForAnnotations -mmq 0 -nda 
INFO  19:44:38,329 HelpFormatter - Executing as sredmond@a2e-bigmem-0.broadinstitute.org on Linux 2.6.32-504.8.1.el6.x86_64 amd64; Java HotSpot(TM) 64-Bit Server VM 1.8.0-b132. 
INFO  19:44:38,330 HelpFormatter - Date/Time: 2015/07/13 19:44:38 
INFO  19:44:38,330 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  19:44:38,331 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  19:44:39,251 GenomeAnalysisEngine - Strictness is SILENT 
INFO  19:44:39,578 GenomeAnalysisEngine - Downsampling Settings: Method: BY_SAMPLE, Target Coverage: 500 
INFO  19:44:39,592 SAMDataSource$SAMReaders - Initializing SAMRecords in serial 
INFO  19:44:39,663 SAMDataSource$SAMReaders - Done initializing BAM readers: total time 0.07 
INFO  19:44:40,542 GATKRunReport - Uploaded run statistics report to AWS S3 
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
Picked up _JAVA_OPTIONS: -Xmx2000m -Xms50m
INFO  19:44:48,236 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  19:44:48,241 HelpFormatter - The Genome Analysis Toolkit (GATK) v3.4-0-g7e26428, Compiled 2015/05/15 03:25:41 
INFO  19:44:48,250 HelpFormatter - Copyright (c) 2010 The Broad Institute 
INFO  19:44:48,250 HelpFormatter - For support and documentation go to http://www.broadinstitute.org/gatk 
INFO  19:44:48,258 HelpFormatter - Program Args: -T HaplotypeCaller -I bams.ITlist -ploidy 1 -R /seq/plasmodium/sredmond//refs/PfIT_v3.fasta -o haplo.rIT.ALL.vcf -gt_mode DISCOVERY -out_mode EMIT_VARIANTS_ONLY -stand_emit_conf 0 -stand_call_conf 30 -L PfIT_07_v3:450000-750000 -APO haplo.rIT.ALL.apo -globalMAPQ -1 -pcrModel NONE -useFilteredReadsForAnnotations -mmq 0 -nda 
INFO  19:44:48,266 HelpFormatter - Executing as sredmond@a2e-bigmem-0.broadinstitute.org on Linux 2.6.32-504.8.1.el6.x86_64 amd64; Java HotSpot(TM) 64-Bit Server VM 1.8.0-b132. 
INFO  19:44:48,266 HelpFormatter - Date/Time: 2015/07/13 19:44:48 
INFO  19:44:48,267 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  19:44:48,267 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  19:44:49,901 GATKRunReport - Uploaded run statistics report to AWS S3 
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
##### ERROR MESSAGE: Invalid command line: The GATK reads argument (-I, --input_file) supports only BAM/CRAM files with the .bam/.cram extension and lists of BAM/CRAM files with the .list extension, but the file bams.ITlist has neither extension.  Please ensure that your BAM/CRAM file or list of BAM/CRAM files is in the correct format, update the extension, and try again.
##### ERROR ------------------------------------------------------------------------------------------
