use strict;
my $windowSize = 20;

open IN,"01_run.pl.regions" or die $!;
my $head = <IN>;
my @regions = ();
my $regionCount = 0;
while (my $line = <IN>)
{
	chomp $line;
	push @regions, $line;
	$regionCount++;
}
close IN;
print "Read $regionCount regions\n";


my @samples = qw (Test);
foreach my $sample (@samples)
{
	my $readCount = 0;
	my $printedCount = 0;
	my $outfile = "$0.$sample.out";
	open OUT, ">$outfile" or die $!;
	my $outCount = 0;
	
	my $currChr = -1;
	my $currStart = -1;
	my $currEnd = -1;
	my @currVals = -1;
	my $currVal = -1;

	print OUT "chr\tstart\tend\texp\texp_noctl\texp_unique\texp_reads\tctl\tctl_unique\tctl_reads\tpct_exp\n";
	foreach my $region (@regions)
	{
		my @lineEls = split "\t", $region;
		my $regChr = $lineEls[0];
		my $regStart = $lineEls[1];
		my $regEnd = $lineEls[2];

		my %expRegionData;
		my %ctlRegionData;
		my %regionKeys;

#		print "opening data/$sample.$chr\_$start\_$end\n";
		my $regionOutput = "data/$sample.$regChr\_$regStart\_$regEnd";
		die "region $regionOutput does not exist" unless -e $regionOutput;
		open IN, $regionOutput or die $!;
		while (my $line = <IN>)
		{
			$readCount++;
			chomp $line;
			my @lineEls = split "\t", $line;
			my $chr = $lineEls[0];
			my $start = $lineEls[1];
			my $end = $start + 2*$windowSize;
			my @vals = @lineEls;
			my $expVal = $lineEls[4]/($lineEls[5]+1);
			my $ctlVal = ($lineEls[7]+1)/($lineEls[8] + 1);
			my $val = $expVal/$ctlVal;

			#if adjacent, just update, don't print
			if ($chr eq $currChr and $start <= $currEnd)
			{
				$currEnd = $end;
				if ($val > $currVal)
				{
					@currVals = @vals;
					$currVal = $val;
				}
			}
			#else print
			else
			{
				if ($currStart > -1)
				{
					$currVals[1] = "$currStart\t$currEnd";
					print OUT join("\t", @currVals)."\t$currVal\n";
					$printedCount++;
				}
				$currChr = $chr;
				$currStart = $start;
				$currEnd = $end;
				@currVals = @vals;
				$currVal = $val;
			}
		}
		$currVals[1] = "$currStart\t$currEnd";
		print OUT join("\t", @currVals)."\t$currVal\n";
		$printedCount++;
	}
	close OUT;
	print "Finished. printed $outfile for $sample\n";
	print "Read $readCount, merged and printed $printedCount\n";
}
print "Finished\n";
