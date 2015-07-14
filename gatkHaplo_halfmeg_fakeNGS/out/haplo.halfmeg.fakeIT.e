Picked up _JAVA_OPTIONS: -Xmx2000m -Xms50m
INFO  20:00:09,919 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  20:00:09,922 HelpFormatter - The Genome Analysis Toolkit (GATK) v3.4-0-g7e26428, Compiled 2015/05/15 03:25:41 
INFO  20:00:09,922 HelpFormatter - Copyright (c) 2010 The Broad Institute 
INFO  20:00:09,923 HelpFormatter - For support and documentation go to http://www.broadinstitute.org/gatk 
INFO  20:00:09,927 HelpFormatter - Program Args: -T HaplotypeCaller -I bams.IT.list -ploidy 1 -R /seq/plasmodium/sredmond//refs/PfIT_v3.fasta -o haplo.rIT.ALL.vcf -gt_mode DISCOVERY -out_mode EMIT_VARIANTS_ONLY -stand_emit_conf 0 -stand_call_conf 30 -L PfIT_07_v3:450000-750000 -APO haplo.rIT.ALL.apo -globalMAPQ -1 -pcrModel NONE -useFilteredReadsForAnnotations -mmq 0 -nda 
INFO  20:00:09,937 HelpFormatter - Executing as sredmond@a2e-bigmem-0.broadinstitute.org on Linux 2.6.32-504.8.1.el6.x86_64 amd64; Java HotSpot(TM) 64-Bit Server VM 1.8.0-b132. 
INFO  20:00:09,937 HelpFormatter - Date/Time: 2015/07/13 20:00:09 
INFO  20:00:09,937 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  20:00:09,938 HelpFormatter - -------------------------------------------------------------------------------- 
INFO  20:00:10,591 GenomeAnalysisEngine - Strictness is SILENT 
INFO  20:00:10,755 GenomeAnalysisEngine - Downsampling Settings: Method: BY_SAMPLE, Target Coverage: 500 
INFO  20:00:10,766 SAMDataSource$SAMReaders - Initializing SAMRecords in serial 
INFO  20:00:10,829 SAMDataSource$SAMReaders - Done initializing BAM readers: total time 0.06 
INFO  20:00:10,854 IntervalUtils - Processing 300001 bp from intervals 
INFO  20:00:10,955 GenomeAnalysisEngine - Preparing for traversal over 3 BAM files 
INFO  20:00:10,994 GenomeAnalysisEngine - Done preparing for traversal 
INFO  20:00:10,994 ProgressMeter - [INITIALIZATION COMPLETE; STARTING PROCESSING] 
INFO  20:00:10,995 ProgressMeter -                 |      processed |    time |         per 1M |           |   total | remaining 
INFO  20:00:10,995 ProgressMeter -        Location | active regions | elapsed | active regions | completed | runtime |   runtime 
INFO  20:00:10,996 HaplotypeCaller - Currently, physical phasing is not available when ploidy is different than 2; therefore it won't be performed 
INFO  20:00:11,102 StrandBiasTest - SAM/BAM data was found. Attempting to use read data to calculate strand bias annotations values. 
INFO  20:00:11,102 StrandBiasTest - SAM/BAM data was found. Attempting to use read data to calculate strand bias annotations values. 
Using SSE4.1 accelerated implementation of PairHMM
INFO  20:00:11,903 VectorLoglessPairHMM - libVectorLoglessPairHMM unpacked successfully from GATK jar file 
INFO  20:00:11,904 VectorLoglessPairHMM - Using vectorized implementation of PairHMM 
WARN  20:00:16,090 HaplotypeScore - Annotation will not be calculated, must be called from UnifiedGenotyper 
WARN  20:00:16,091 InbreedingCoeff - Annotation will not be calculated, must provide a valid PED file (-ped) from the command line. 
WARN  20:00:16,094 AnnotationUtils - Annotation will not be calculated, genotype is not called 
WARN  20:00:35,601 ExactAFCalculator - this tool is currently set to genotype at most 6 alternate alleles in a given context, but the context at PfIT_07_v3:484281 has 8 alternate alleles so only the top alleles will be used; see the --max_alternate_alleles argument 
INFO  20:00:41,000 ProgressMeter - PfIT_07_v3:485616              0.0    30.0 s           49.6 w       11.9%     4.2 m       3.7 m 
INFO  20:01:11,001 ProgressMeter - PfIT_07_v3:491675              0.0    60.0 s           99.2 w       13.9%     7.2 m       6.2 m 
INFO  20:01:41,003 ProgressMeter - PfIT_07_v3:494038              0.0    90.0 s          148.8 w       14.7%    10.2 m       8.7 m 
INFO  20:02:11,004 ProgressMeter - PfIT_07_v3:496361              0.0   120.0 s          198.4 w       15.5%    12.9 m      10.9 m 
INFO  20:02:41,006 ProgressMeter - PfIT_07_v3:499077              0.0     2.5 m          248.0 w       16.4%    15.3 m      12.8 m 
WARN  20:02:53,717 ExactAFCalculator - this tool is currently set to genotype at most 6 alternate alleles in a given context, but the context at PfIT_07_v3:499675 has 7 alternate alleles so only the top alleles will be used; see the --max_alternate_alleles argument 
INFO  20:03:11,007 ProgressMeter - PfIT_07_v3:503135              0.0     3.0 m          297.6 w       17.7%    16.9 m      13.9 m 
INFO  20:03:41,009 ProgressMeter - PfIT_07_v3:505192              0.0     3.5 m          347.2 w       18.4%    19.0 m      15.5 m 
INFO  20:04:11,010 ProgressMeter - PfIT_07_v3:507735              0.0     4.0 m          396.9 w       19.2%    20.8 m      16.8 m 
INFO  20:04:41,012 ProgressMeter - PfIT_07_v3:509799              0.0     4.5 m          446.5 w       19.9%    22.6 m      18.1 m 
INFO  20:05:11,014 ProgressMeter - PfIT_07_v3:511385              0.0     5.0 m          496.1 w       20.5%    24.4 m      19.4 m 
INFO  20:05:41,015 ProgressMeter - PfIT_07_v3:518618              0.0     5.5 m          545.7 w       22.9%    24.0 m      18.5 m 
INFO  20:06:11,017 ProgressMeter - PfIT_07_v3:522411              0.0     6.0 m          595.3 w       24.1%    24.9 m      18.9 m 
INFO  20:06:41,018 ProgressMeter - PfIT_07_v3:524656              0.0     6.5 m          644.9 w       24.9%    26.1 m      19.6 m 
INFO  20:07:11,021 ProgressMeter - PfIT_07_v3:526291              0.0     7.0 m          694.5 w       25.4%    27.5 m      20.5 m 
INFO  20:07:41,023 ProgressMeter - PfIT_07_v3:529204              0.0     7.5 m          744.1 w       26.4%    28.4 m      20.9 m 
WARN  20:08:07,726 ExactAFCalculator - this tool is currently set to genotype at most 6 alternate alleles in a given context, but the context at PfIT_07_v3:532958 has 9 alternate alleles so only the top alleles will be used; see the --max_alternate_alleles argument 
INFO  20:08:11,024 ProgressMeter - PfIT_07_v3:533619              0.0     8.0 m          793.7 w       27.9%    28.7 m      20.7 m 
INFO  20:08:41,026 ProgressMeter - PfIT_07_v3:534988              0.0     8.5 m          843.3 w       28.3%    30.0 m      21.5 m 
INFO  20:09:11,027 ProgressMeter - PfIT_07_v3:538227              0.0     9.0 m          892.9 w       29.4%    30.6 m      21.6 m 
INFO  20:09:41,029 ProgressMeter - PfIT_07_v3:539437              0.0     9.5 m          942.5 w       29.8%    31.9 m      22.4 m 
INFO  20:10:11,030 ProgressMeter - PfIT_07_v3:541595              0.0    10.0 m          992.1 w       30.5%    32.8 m      22.8 m 
INFO  20:10:41,032 ProgressMeter - PfIT_07_v3:542764              0.0    10.5 m         1041.7 w       30.9%    34.0 m      23.5 m 
INFO  20:11:11,034 ProgressMeter - PfIT_07_v3:546048              0.0    11.0 m         1091.3 w       32.0%    34.4 m      23.4 m 
INFO  20:11:41,036 ProgressMeter - PfIT_07_v3:552111              0.0    11.5 m         1140.9 w       34.0%    33.8 m      22.3 m 
INFO  20:12:11,037 ProgressMeter - PfIT_07_v3:553060              0.0    12.0 m         1190.5 w       34.4%    34.9 m      22.9 m 
INFO  20:12:41,039 ProgressMeter - PfIT_07_v3:553776              0.0    12.5 m         1240.2 w       34.6%    36.1 m      23.6 m 
INFO  20:13:11,041 ProgressMeter - PfIT_07_v3:557102              0.0    13.0 m         1289.8 w       35.7%    36.4 m      23.4 m 
INFO  20:13:41,043 ProgressMeter - PfIT_07_v3:564011              0.0    13.5 m         1339.4 w       38.0%    35.5 m      22.0 m 
INFO  20:14:11,045 ProgressMeter - PfIT_07_v3:574838              0.0    14.0 m         1389.0 w       41.6%    33.6 m      19.6 m 
INFO  20:14:41,047 ProgressMeter - PfIT_07_v3:578067              0.0    14.5 m         1438.6 w       42.7%    34.0 m      19.5 m 
INFO  20:15:11,048 ProgressMeter - PfIT_07_v3:581467              0.0    15.0 m         1488.2 w       43.8%    34.2 m      19.2 m 
INFO  20:15:41,050 ProgressMeter - PfIT_07_v3:582841              0.0    15.5 m         1537.8 w       44.3%    35.0 m      19.5 m 
INFO  20:16:11,052 ProgressMeter - PfIT_07_v3:585636              0.0    16.0 m         1587.4 w       45.2%    35.4 m      19.4 m 
INFO  20:16:41,056 ProgressMeter - PfIT_07_v3:586764              0.0    16.5 m         1637.0 w       45.6%    36.2 m      19.7 m 
INFO  20:17:11,058 ProgressMeter - PfIT_07_v3:588136              0.0    17.0 m         1686.6 w       46.0%    36.9 m      19.9 m 
INFO  20:17:41,060 ProgressMeter - PfIT_07_v3:588945              0.0    17.5 m         1736.2 w       46.3%    37.8 m      20.3 m 
INFO  20:18:11,062 ProgressMeter - PfIT_07_v3:590460              0.0    18.0 m         1785.8 w       46.8%    38.4 m      20.4 m 
INFO  20:18:41,064 ProgressMeter - PfIT_07_v3:591171              0.0    18.5 m         1835.4 w       47.1%    39.3 m      20.8 m 
INFO  20:19:11,065 ProgressMeter - PfIT_07_v3:592198              0.0    19.0 m         1885.0 w       47.4%    40.1 m      21.1 m 
INFO  20:19:41,067 ProgressMeter - PfIT_07_v3:593625              0.0    19.5 m         1934.6 w       47.9%    40.7 m      21.2 m 
INFO  20:20:11,068 ProgressMeter - PfIT_07_v3:595985              0.0    20.0 m         1984.2 w       48.7%    41.1 m      21.1 m 
INFO  20:20:41,070 ProgressMeter - PfIT_07_v3:598154              0.0    20.5 m         2033.9 w       49.4%    41.5 m      21.0 m 
INFO  20:21:11,071 ProgressMeter - PfIT_07_v3:600643              0.0    21.0 m         2083.5 w       50.2%    41.8 m      20.8 m 
INFO  20:21:41,073 ProgressMeter - PfIT_07_v3:604102              0.0    21.5 m         2133.1 w       51.4%    41.9 m      20.4 m 
INFO  20:22:11,074 ProgressMeter - PfIT_07_v3:607336              0.0    22.0 m         2182.7 w       52.4%    41.9 m      19.9 m 
INFO  20:22:41,076 ProgressMeter - PfIT_07_v3:611391              0.0    22.5 m         2232.3 w       53.8%    41.8 m      19.3 m 
INFO  20:23:11,077 ProgressMeter - PfIT_07_v3:613718              0.0    23.0 m         2281.9 w       54.6%    42.1 m      19.1 m 
INFO  20:23:41,088 ProgressMeter - PfIT_07_v3:616974              0.0    23.5 m         2331.5 w       55.7%    42.2 m      18.7 m 
INFO  20:24:11,090 ProgressMeter - PfIT_07_v3:621179              0.0    24.0 m         2381.1 w       57.1%    42.1 m      18.1 m 
INFO  20:24:41,091 ProgressMeter - PfIT_07_v3:625348              0.0    24.5 m         2430.7 w       58.4%    41.9 m      17.4 m 
INFO  20:25:11,092 ProgressMeter - PfIT_07_v3:628627              0.0    25.0 m         2480.3 w       59.5%    42.0 m      17.0 m 
INFO  20:25:41,094 ProgressMeter - PfIT_07_v3:631862              0.0    25.5 m         2529.9 w       60.6%    42.1 m      16.6 m 
INFO  20:26:11,096 ProgressMeter - PfIT_07_v3:635794              0.0    26.0 m         2579.5 w       61.9%    42.0 m      16.0 m 
INFO  20:26:41,098 ProgressMeter - PfIT_07_v3:637461              0.0    26.5 m         2629.1 w       62.5%    42.4 m      15.9 m 
INFO  20:27:11,099 ProgressMeter - PfIT_07_v3:642858              0.0    27.0 m         2678.7 w       64.3%    42.0 m      15.0 m 
INFO  20:27:41,100 ProgressMeter - PfIT_07_v3:650801              0.0    27.5 m         2728.3 w       66.9%    41.1 m      13.6 m 
INFO  20:28:11,101 ProgressMeter - PfIT_07_v3:653827              0.0    28.0 m         2778.0 w       67.9%    41.2 m      13.2 m 
INFO  20:28:41,103 ProgressMeter - PfIT_07_v3:654946              0.0    28.5 m         2827.6 w       68.3%    41.7 m      13.2 m 
INFO  20:29:11,104 ProgressMeter - PfIT_07_v3:658150              0.0    29.0 m         2877.2 w       69.4%    41.8 m      12.8 m 
INFO  20:29:41,138 ProgressMeter - PfIT_07_v3:660555              0.0    29.5 m         2926.8 w       70.2%    42.0 m      12.5 m 
INFO  20:30:11,139 ProgressMeter - PfIT_07_v3:663135              0.0    30.0 m         2976.4 w       71.0%    42.2 m      12.2 m 
INFO  20:30:41,140 ProgressMeter - PfIT_07_v3:667046              0.0    30.5 m         3026.0 w       72.3%    42.2 m      11.7 m 
INFO  20:31:11,142 ProgressMeter - PfIT_07_v3:669824              0.0    31.0 m         3075.6 w       73.3%    42.3 m      11.3 m 
INFO  20:31:41,143 ProgressMeter - PfIT_07_v3:672142              0.0    31.5 m         3125.2 w       74.0%    42.5 m      11.0 m 
INFO  20:32:11,144 ProgressMeter - PfIT_07_v3:674884              0.0    32.0 m         3174.9 w       75.0%    42.7 m      10.7 m 
INFO  20:32:41,146 ProgressMeter - PfIT_07_v3:677378              0.0    32.5 m         3224.5 w       75.8%    42.9 m      10.4 m 
INFO  20:33:11,147 ProgressMeter - PfIT_07_v3:680937              0.0    33.0 m         3274.1 w       77.0%    42.9 m       9.9 m 
INFO  20:33:41,149 ProgressMeter - PfIT_07_v3:683930              0.0    33.5 m         3323.7 w       78.0%    43.0 m       9.5 m 
INFO  20:34:11,150 ProgressMeter - PfIT_07_v3:685142              0.0    34.0 m         3373.3 w       78.4%    43.4 m       9.4 m 
INFO  20:34:41,151 ProgressMeter - PfIT_07_v3:688335              0.0    34.5 m         3422.9 w       79.4%    43.4 m       8.9 m 
INFO  20:35:11,153 ProgressMeter - PfIT_07_v3:690394              0.0    35.0 m         3472.5 w       80.1%    43.7 m       8.7 m 
INFO  20:35:41,154 ProgressMeter - PfIT_07_v3:693464              0.0    35.5 m         3522.1 w       81.2%    43.7 m       8.2 m 
INFO  20:36:11,156 ProgressMeter - PfIT_07_v3:694234              0.0    36.0 m         3571.7 w       81.4%    44.2 m       8.2 m 
WARN  20:36:32,932 ExactAFCalculator - this tool is currently set to genotype at most 6 alternate alleles in a given context, but the context at PfIT_07_v3:694817 has 14 alternate alleles so only the top alleles will be used; see the --max_alternate_alleles argument 
INFO  20:36:41,158 ProgressMeter - PfIT_07_v3:697675              0.0    36.5 m         3621.3 w       82.6%    44.2 m       7.7 m 
INFO  20:37:11,159 ProgressMeter - PfIT_07_v3:700505              0.0    37.0 m         3670.9 w       83.5%    44.3 m       7.3 m 
INFO  20:37:41,161 ProgressMeter - PfIT_07_v3:703103              0.0    37.5 m         3720.5 w       84.4%    44.4 m       6.9 m 
WARN  20:38:06,153 ExactAFCalculator - this tool is currently set to genotype at most 6 alternate alleles in a given context, but the context at PfIT_07_v3:703433 has 7 alternate alleles so only the top alleles will be used; see the --max_alternate_alleles argument 
INFO  20:38:11,162 ProgressMeter - PfIT_07_v3:704406              0.0    38.0 m         3770.1 w       84.8%    44.8 m       6.8 m 
WARN  20:38:18,128 ExactAFCalculator - this tool is currently set to genotype at most 6 alternate alleles in a given context, but the context at PfIT_07_v3:704055 has 8 alternate alleles so only the top alleles will be used; see the --max_alternate_alleles argument 
INFO  20:38:41,163 ProgressMeter - PfIT_07_v3:706041              0.0    38.5 m         3819.7 w       85.3%    45.1 m       6.6 m 
WARN  20:38:58,430 ExactAFCalculator - this tool is currently set to genotype at most 6 alternate alleles in a given context, but the context at PfIT_07_v3:706781 has 7 alternate alleles so only the top alleles will be used; see the --max_alternate_alleles argument 
INFO  20:39:11,165 ProgressMeter - PfIT_07_v3:708155              0.0    39.0 m         3869.3 w       86.1%    45.3 m       6.3 m 
INFO  20:39:41,166 ProgressMeter - PfIT_07_v3:710168              0.0    39.5 m         3918.9 w       86.7%    45.5 m       6.0 m 
WARN  20:39:48,921 ExactAFCalculator - this tool is currently set to genotype at most 6 alternate alleles in a given context, but the context at PfIT_07_v3:710167 has 7 alternate alleles so only the top alleles will be used; see the --max_alternate_alleles argument 
INFO  20:40:11,167 ProgressMeter - PfIT_07_v3:711332              0.0    40.0 m         3968.5 w       87.1%    45.9 m       5.9 m 
INFO  20:40:41,169 ProgressMeter - PfIT_07_v3:713018              0.0    40.5 m         4018.1 w       87.7%    46.2 m       5.7 m 
INFO  20:41:11,170 ProgressMeter - PfIT_07_v3:714483              0.0    41.0 m         4067.8 w       88.2%    46.5 m       5.5 m 
INFO  20:41:41,171 ProgressMeter - PfIT_07_v3:717853              0.0    41.5 m         4117.4 w       89.3%    46.5 m       5.0 m 
INFO  20:42:11,173 ProgressMeter - PfIT_07_v3:719538              0.0    42.0 m         4167.0 w       89.8%    46.7 m       4.7 m 
INFO  20:42:41,174 ProgressMeter - PfIT_07_v3:721197              0.0    42.5 m         4216.6 w       90.4%    47.0 m       4.5 m 
INFO  20:43:11,195 ProgressMeter - PfIT_07_v3:723380              0.0    43.0 m         4266.2 w       91.1%    47.2 m       4.2 m 
INFO  20:43:41,197 ProgressMeter - PfIT_07_v3:725613              0.0    43.5 m         4315.8 w       91.9%    47.3 m       3.8 m 
WARN  20:43:44,404 ExactAFCalculator - this tool is currently set to genotype at most 6 alternate alleles in a given context, but the context at PfIT_07_v3:725297 has 7 alternate alleles so only the top alleles will be used; see the --max_alternate_alleles argument 
INFO  20:44:11,198 ProgressMeter - PfIT_07_v3:726285              0.0    44.0 m         4365.4 w       92.1%    47.8 m       3.8 m 
WARN  20:44:20,600 ExactAFCalculator - this tool is currently set to genotype at most 6 alternate alleles in a given context, but the context at PfIT_07_v3:726380 has 8 alternate alleles so only the top alleles will be used; see the --max_alternate_alleles argument 
INFO  20:44:41,200 ProgressMeter - PfIT_07_v3:727604              0.0    44.5 m         4415.0 w       92.5%    48.1 m       3.6 m 
INFO  20:45:11,201 ProgressMeter - PfIT_07_v3:729206              0.0    45.0 m         4464.6 w       93.1%    48.4 m       3.4 m 
INFO  20:45:41,202 ProgressMeter - PfIT_07_v3:731753              0.0    45.5 m         4514.2 w       93.9%    48.4 m       2.9 m 
WARN  20:45:49,578 ExactAFCalculator - this tool is currently set to genotype at most 6 alternate alleles in a given context, but the context at PfIT_07_v3:731454 has 10 alternate alleles so only the top alleles will be used; see the --max_alternate_alleles argument 
INFO  20:46:11,204 ProgressMeter - PfIT_07_v3:732229              0.0    46.0 m         4563.8 w       94.1%    48.9 m       2.9 m 
WARN  20:46:20,577 ExactAFCalculator - this tool is currently set to genotype at most 6 alternate alleles in a given context, but the context at PfIT_07_v3:731947 has 12 alternate alleles so only the top alleles will be used; see the --max_alternate_alleles argument 
WARN  20:46:20,582 ExactAFCalculator - this tool is currently set to genotype at most 6 alternate alleles in a given context, but the context at PfIT_07_v3:731991 has 10 alternate alleles so only the top alleles will be used; see the --max_alternate_alleles argument 
INFO  20:46:41,205 ProgressMeter - PfIT_07_v3:733435              0.0    46.5 m         4613.4 w       94.5%    49.2 m       2.7 m 
INFO  20:47:11,207 ProgressMeter - PfIT_07_v3:734565              0.0    47.0 m         4663.0 w       94.9%    49.5 m       2.5 m 
INFO  20:47:41,208 ProgressMeter - PfIT_07_v3:737574              0.0    47.5 m         4712.7 w       95.9%    49.6 m       2.1 m 
INFO  20:48:11,217 ProgressMeter - PfIT_07_v3:741516              0.0    48.0 m         4762.3 w       97.2%    49.4 m      83.0 s 
INFO  20:48:41,219 ProgressMeter - PfIT_07_v3:745640              0.0    48.5 m         4811.9 w       98.5%    49.2 m      42.0 s 
INFO  20:48:58,741 VectorLoglessPairHMM - Time spent in setup for JNI call : 0.393171741 
INFO  20:48:58,742 PairHMM - Total compute time in PairHMM computeLikelihoods() : 2397.237217112 
INFO  20:48:58,742 HaplotypeCaller - Ran local assembly on 625 active regions 
INFO  20:48:58,774 ProgressMeter -            done         300001.0    48.8 m            2.7 h      100.0%    48.8 m       0.0 s 
INFO  20:48:58,775 ProgressMeter - Total runtime 2927.78 secs, 48.80 min, 0.81 hours 
INFO  20:48:58,775 MicroScheduler - 1789 reads were filtered out during the traversal out of approximately 300820 total reads (0.59%) 
INFO  20:48:58,776 MicroScheduler -   -> 0 reads (0.00% of total) failing BadCigarFilter 
INFO  20:48:58,777 MicroScheduler -   -> 0 reads (0.00% of total) failing DuplicateReadFilter 
INFO  20:48:58,777 MicroScheduler -   -> 0 reads (0.00% of total) failing FailsVendorQualityCheckFilter 
INFO  20:48:58,777 MicroScheduler -   -> 0 reads (0.00% of total) failing HCMappingQualityFilter 
INFO  20:48:58,778 MicroScheduler -   -> 0 reads (0.00% of total) failing MalformedReadFilter 
INFO  20:48:58,778 MicroScheduler -   -> 0 reads (0.00% of total) failing MappingQualityUnavailableFilter 
INFO  20:48:58,778 MicroScheduler -   -> 1702 reads (0.57% of total) failing NotPrimaryAlignmentFilter 
INFO  20:48:58,778 MicroScheduler -   -> 87 reads (0.03% of total) failing UnmappedReadFilter 
INFO  20:48:59,603 GATKRunReport - Uploaded run statistics report to AWS S3 
