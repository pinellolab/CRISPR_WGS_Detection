use strict;
my $usage = "Usage: perl $0 {bp to pad} {casoffinder output} {output file} {cleavage offset}";
#number of basepairs to pad
my $padding = shift || 42;
my $casoffinderFile = shift;
my $outputFile = shift || "$casoffinderFile.pad$padding.txt";
my $cleavageOffset = shift || -3; #-3 for cas9
die "Cannot find casoffinder file $casoffinderFile\n$usage" unless -e $casoffinderFile;
print "Running with parameters:
padding: $padding
casoffinderFile: $casoffinderFile
outputFile: $outputFile
cleavageOffset: $cleavageOffset\n";


open IN, $casoffinderFile or die $!;
my $head = <IN>;
my $bulgeMode = 0;
if ($head =~ /^#Bulge type/)
{
	$bulgeMode = 1;
}
seek(IN,0,0);

open OUT, ">$outputFile" or die $!;
my $pamoutfile = "$outputFile.pam.txt";
open PAM, ">$pamoutfile" or die $!;
my $count = 0;
while (my $line = <IN>)
{
	chomp $line;
	my @lineEls = split "\t", $line;
	my $chr = $lineEls[1];
	my $start = $lineEls[2];
	my $guideWithPAM = $lineEls[3];
	my $mismatch = $lineEls[5];
	my $strand = $lineEls[4];
	my $newName = "offby$mismatch"."_$count";
    my $guideseq = $lineEls[9];

	if ($bulgeMode)
	{
		$chr = $lineEls[3];
		$start = $lineEls[4];
		$guideWithPAM = $lineEls[2];
		$mismatch = $lineEls[6];
		$strand = $lineEls[5];
		$newName = "$lineEls[0]$lineEls[7]$lineEls[6]"."_$count";
        $guideseq = $lineEls[9];
	}

	my $guide = substr($guideWithPAM,0,-3);
	$guide =~ s/-//g;
	my $cleavagePos = $start + length($guide) + $cleavageOffset;
	if ($strand eq "-")
	{
		$cleavagePos = $start + (-1*$cleavageOffset) + 3; #+3 for PAM
	}

	my $newStart = $cleavagePos - $padding;
	my $newEnd = $cleavagePos + $padding + 1;

	print OUT "$chr\t$newStart\t$newEnd\t$newName\t$guide\t$guideseq\n";
	print PAM "$chr\t$newStart\t$newEnd\t$newName\t$guide\t$guideWithPAM\t$line\t$guideseq\n";
	$count++;

}
print "Finished. Printed $count to $outputFile\n";
