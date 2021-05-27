use strict;
my @samples = qw (Test);
mkdir 'crispresso' unless -d 'crispresso';
mkdir 'crispressoControl' unless -d 'crispressoControl';
my $numToRun = 20;

my $expReadsCol = 6;
my $ctlReadsCol = 9;
my $pctCol = 10;
my $ctlBam = "../../../ctl.bam";

open CMD, ">$0.run.sh" or die $!;
foreach my $sample (@samples)
{
	open IN, "02_aggregate.pl.$sample.out" or die $!;
	my %data = ();
	my $head = <IN>;
	chomp $head;
	my @headEls = split "\t", $head;
	if ($headEls[$expReadsCol] ne "exp_reads")
	{
		die "expecting exp_reads at $headEls[$expReadsCol]\n";
	}
	if ($headEls[$ctlReadsCol] ne "ctl_reads")
	{
		die "expecting ctl_reads at $headEls[$ctlReadsCol]\n";
	}
	if ($headEls[$pctCol] ne "pct_exp")
	{
		die "expecting pct_exp at $headEls[$pctCol]\n";
	}
	my $readLinesCount = 0;
	my $keptLinesCount = 0;
	while (my $line = <IN>)
	{
		$readLinesCount++;
		chomp $line;
		my @lineEls = split "\t", $line;
		my $expReads = $lineEls[$expReadsCol];
		my $ctlReads = $lineEls[$ctlReadsCol];
		my $pct = $lineEls[$pctCol];

		if ($expReads > 100 and $ctlReads > 100 and $pct > 3)
		{
			$data{$line} = $pct;
			$keptLinesCount++;
		}
	}
	print "for $sample kept $keptLinesCount/$readLinesCount\n";


	my $outfile = "$0.$sample.regions";
	open OUT, ">$outfile" or die $!;
	my @sortedVals = sort {$data{$b} <=> $data{$a} } keys %data;
	for (my $i = 0; $i < $numToRun; $i++)
	{
		my $val = $sortedVals[$i];
		chomp $val;
		my @valEls = split "\t", $val;

		my $start = $valEls[1];
		my $end = $start + 40;
		my $name = "hit_$i"."_".int($data{$val}+.5);
		print OUT "$valEls[0]\t$valEls[1]\t$valEls[2]\t$name\n";
	}
	my $params = " --exclude_bp_from_right 0 --exclude_bp_from_left 0 --plot_window_size 10 -p $numToRun --no_rerun --default_min_aln_score 10 ";
	print CMD "CRISPRessoWGS -b ../../../$sample.bam -f $0.$sample.regions -o crispresso -r GENOMES/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa $params\n";
	print CMD "CRISPRessoWGS -b $ctlBam -f $0.$sample.regions -o crispressoControl -r GENOMES/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa $params\n";
}


	
print "Finished\n";
print "run $0.run.sh\n";

