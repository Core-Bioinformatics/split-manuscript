#########################################################
#PREPROCESSING
#########################################################
# 1.1 transform redundant fasta file into non-redundant fasta file
#if a non-redundant file is given as input the script will return the same file to output

use Bio::SeqIO;

$usage = "perl R2NR.pl fastaFile out";
$fasta = shift or die $usage;
$out   = shift or die $usage;

print "R2NR:: input file: $fasta \n";

my $fa = Bio::SeqIO->new('-format' => 'fasta', '-file' =>  "$fasta");

open out, '>'.$out
	or die "cannot open output";

$s = '';
$count = 0;
$cnt = 0;
$initial = 1;
while(my $seq = $fa->next_seq()) 
{
	my $sequence = $seq->seq();        # get the sequence as a string
	my $id       = $seq->display_id();  # identifier
#	@abun = split(/\(/,$id);
#	@abn  = split(/\)/,$abun[1]);
	
	push @seqs, [$id, $sequence];
}

print "total number of sequences read : $#seqs \n";
@ss = sort{$a->[1] cmp $b->[1]} @seqs;

$s = '';
$count = 0;
$cnt = 0;
$cntNR = 0;
$initial = 1;

for($i = 0; $i <= $#ss; $i++)
{
	if($s ne $ss[$i][1])
	{
		if($initial == 1)
		{
			$initial = 0;
		}
		else
		{
			print out ">$s-$count \n";
			print out "$s\n";
			$cntNR++;
		}
		$cnt += $count;
		$s = $ss[$i][1];
		
#		$defined = 0;
#		$defined = 1 if defined $abn[0];
#		if($defined)
#		{
#			$count = $abn[0];
#		}
#		else
		{
			$count = 1;
		}
	}
	else
	{
#		$defined = 0;
#		$defined = 1 if defined $abn[0];
#		if($defined)
#		{
#			$count += $abn[0];
#		}
#		else
		{
			$count ++;
		}
	}
}

print "R2NR:: total redundant sequences: $cnt \n";
print "R2NR:: total non-redundant sequences: $cntNR \n";

exit;