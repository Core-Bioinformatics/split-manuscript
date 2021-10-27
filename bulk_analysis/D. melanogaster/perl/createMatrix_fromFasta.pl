#creates matrix of exprssion from non redundant fa files
use Bio::SeqIO;
$usage = "perl createMatrixExpression.pl inFiles out";

$inF = shift or die $usage;
$out = shift or die $usage;

print "all files : $inF \n";
print "out : $out \n";

open inF, $inF
	or die "cannt open inFiles";
$index = 0;
while(<inF>)
{
	chomp;
	$files[$index][0] = $_;
	#@sp = split/\./;
	#$files[$index][1] = $sp[0];
	$files[$index][1] = $_;

	$index++;
}
print "working with $index files\n";

for($i = 0; $i < $index; $i++)
{
	$fasta = $files[$i][1];
	print "input file : $fasta \n";

	my $fa = Bio::SeqIO->new('-format' => 'fasta', '-file' =>  "$fasta");

	while(my $seq = $fa->next_seq()) 
	{
		my $sequence = $seq->seq();
		my $id       = $seq->display_id();
		@idd = split(/-/,$id);
		#@idd = split(/\(/,$id);
		#@cnt = split(/\)/,$idd[1]);
		push @all,[($i, $sequence, $idd[1])];
	}
}

print "input read \n";

@allS = sort {lc($a->[1]) cmp lc($b->[1])} @all;
print "sorted $#allS \n";

#create series;
$seq = $allS[0][1];
$indS = 0;
for($j = 0; $j < $index; $j++)
{
	$series[$indS][$j] = 0;
}

for($i = 0; $i <= $#allS; $i++)
{
	if($seq eq $allS[$i][1])
	{
		#$series[$indS][0] = $seq;
		$series[$indS][$allS[$i][0]+1] = $allS[$i][2];
	}
	else
	{
		$indS++;
		for($j = 0; $j < $index; $j++)
		{
			$series[$indS][$j] = 0;
		}
		$seq = $allS[$i][1];
		$series[$indS][0] = $seq;
		$series[$indS][$allS[$i][0]+1] = $allS[$i][2];
		
	}
}

print "total series : $indS \n";
print "points in the series : $index \n";

open out, '>'.$out
	or die "cannot open output file";
for($i = 0; $i < $indS; $i++)
{
	for($j = 0; $j <$index; $j++)
	{
		print out "$series[$i][$j], ";
	}
	print out "$series[$i][$j]\n";
}


exit;