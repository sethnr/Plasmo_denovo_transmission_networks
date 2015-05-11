#####
# Job submission wrappers:
#####

#def submit single-job job
sub _submitJob {
    my $command = shift;
    my $dependency = shift;
    
    $dependency = "-w \"".$dependency."\"" if $dependency;
    #print $dependency."\n" if $dependency;
    $dependency = "" if !$dependency;
    my $preexec_com = "";
    $preexec_com = "-E ".$preexec if $preexec;

  # bsub, write to log file...
    my $job_name = "pbj_".$dataset_name."_".$stage;
    #print $job_name."\n";
    my $command = <<"END";
 bsub -J $job_name -o $job_name.%J.out -e $job_name.%J.err \\
 -q $queue -R \'select[mem>$LOCAL_MEM] rusage[mem=$LOCAL_MEM] span[ptile=4]\' \\
 $dependency $preexec_com  -M$LOCAL_MEM -n4  \\
 $jelly $stage $tempfile_sm
END
    $command .= " -x $extras" if $extras;
    _bsub($command);
}


#def submit multi-job job
sub _submitDistJob {
    my $stage = shift;
    my $queue = shift;
    my $dependency = shift;
    
    $dependency = "-w \"".$dependency."\"" if $dependency;
    #print $dependency."\n" if $dependency;
    $dependency = "" if !$dependency;
    my $preexec_com = "";
    $preexec_com = "-E ".$preexec if ($preexec ne "");

  # bsub, write to log file...
    my $job_name = "pbj_".$dataset_name."_".$stage;
    #print $job_name."\n";
    my $command = <<"END";
 bsub -J $job_name -o $job_name.%J.out -e $job_name.%J.err \\
 -q $queue -R \'select[mem>$DIST_MEM_HEAD] rusage[mem=$DIST_MEM_HEAD] span[ptile=4]\' \\
 $dependency $preexec_com  -M$DIST_MEM_HEAD -n4  \\
 $jelly $stage $tempfile_lg
END

    $command .= " -x $extras" if $extras;
    _bsub($command);

    #wait for job to return (all child jobs submitted) 
    my $result = _wait_done();
}

#bsub job and wait till job finishes (collect child jobs after) 
sub _wait_done {
    my $bsub = "";
    my $status = "";
    #count how many wait cycles...
    my $waits = 0;
    print STDERR "  job ".$LASTID." waiting";
    do {
	print STDERR ".";
	$waits++;
	sleep(30);
	$bsub = `bjobs -noheader $LASTID`;
	if($bsub) {
	    my @F = split(/\s+/,$bsub);
	    $status = $F[2];
	}
	#if still waiting and nothing on job list - assume failed to submit
	if ($waits > 3 && $status ne "PEND" 
	    && $status ne "RUN" && $status ne "DONE") {
	    print STDERR "job $LASTID not returning status [$status] - not submitted?\n";
	    exit(1);
	}
	# check through log file for failed jobs;
	# will fail if finds one with EXIT status
	_get_all_jobs(); 
	
    } while ($status ne "EXIT" && $status ne "DONE");
    print STDERR "\n";
    return $status;
}

sub _bsub {
    my $command = shift;
    my $output = `$command`;
    $command =~ s/\n/ /gi;
    if($output =~ m/\<(\d+)\>.*submitted*/gi) {
	$LASTID = $1;
	print LOGFILE $LASTID."\t".$LASTID."\t".$command."\n";
	return $output;
    }
    else{
	print LOGFILE "-1\t-1\tcould not submit command\n".$command;
	$LASTID = "-1";}
}

#make dependency line to wait for end of all current jobs
sub _make_depends_children {
    my @depends = ();
    foreach my $child (@_) {
	my $depend = "done(".$child.")";
	push(@depends,$depend);
    }
    return join(" && ",@depends);
}

######
# get LSF info
######

sub _get_child_jobs {
    my $parent_id = shift;
    my @children;
    open(READLOG,$LOGFILE);
    while(<READLOG>) {
#	print $_;
	my @F = split("\t",$_);
	my ($parent, $child) = @F[0..1];
	if ($parent == -1) {print STDERR "job not submitted, exiting\n";
			    exit(1);
	}
	if ($parent == $parent_id) {
	    push(@children, $child);
	}
    }
    close(READLOG);
    return \@children;
}

sub _get_all_jobs {
    my @jobs;
    open(READLOG,$LOGFILE);
    while(<READLOG>) {
	my @F = split("\t",$_);
	my ($parent, $job) = @F[0..1];
	if ($job == -1) {print STDERR "job not submitted, exiting\n";
			 exit(1);}
	my $status = _get_job_status($job);
	if ($status eq "EXIT") {
	    print STDERR "job ".$job." finished with status ".$status.", exiting\n";
	    system('cat '.$run_folder.'/*'.$job.'*.err');
	    exit(1);}
#	print join("\t","ALLJOBS",$parent, $job, $status)."\n";
	push(@jobs, $job);
    }
    close(READLOG);
    return \@jobs;
}

sub _get_job_status {
    my $job_id = shift;
    my $bjob = `bjobs -noheader $job_id`;
    my $status = "NULL";
    chomp($bjob);
    if(length($bjob) > 0) {
	my @F = split(/\s+/,$bjob);
	$status = $F[2];
    }
    return $status;
}


########
# the stuff that does stuff
########

my $children;

#SETUP: make suffixarray...
  print STDERR "submitting SETUP\n";
  _submitLocalJob("setup","normal");
  print STDERR "  job ".$LASTID." submitted\n";

#MAPPING
  print STDERR "submitting MAPPING (dist)\n";
  _submitDistJob("mapping","normal","done($LASTID)");
  $children = _get_child_jobs($LASTID);
  $LASTDEPENDENCY = _make_depends_children(@{$children});
  print STDERR "  job ".$LASTID." submitted\n";


#SUPPORT
  print STDERR "submitting SUPPORT\n";
  _submitLocalJob("support","normal",$LASTDEPENDENCY);
  print STDERR "  job ".$LASTID." submitted\n";


#EXTRACT (single)
  print STDERR "submitting EXTRACT\n";
  _submitLocalJob("extraction","normal","done($LASTID)");
  print STDERR "  job ".$LASTID." submitted\n";

#ASSEMBLY (multi)
  print STDERR "submitting ASSEMBLY (dist)\n";
  _submitDistJob("assembly","normal","done($LASTID)");
  $children = _get_child_jobs($LASTID);
  $LASTDEPENDENCY = _make_depends_children(@{$children});
  print STDERR "  job ".$LASTID." submitted\n";

#OUTPUT (single)
  print STDERR "submitting OUTPUT\n";
  _submitLocalJob("output","normal",$LASTDEPENDENCY);
  print STDERR "  job ".$LASTID." submitted\n";

sleep(60);
my $all_jobs = _get_all_jobs();

print STDERR "all jobs submitted:\n";
my $get_all = "bjobs ".join(" ",@{$all_jobs});
