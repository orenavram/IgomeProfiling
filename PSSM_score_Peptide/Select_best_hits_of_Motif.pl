use strict;
if (@ARGV<2){die "USAGE: perl $0 <Motifs_Hit_File> <Out>\n";}
my $Motifs_Hits_File=$ARGV[0];
my $Out_File=$ARGV[1];

my @motifs=();
my %Motifs_Name_To_Header=();
my %Motif_Best_Hits=(); # key1: peptide seq; key2: {BEST_MOTIVE | BEST_SCORE }
my %Motivs_hits=(); # key1: motive_name; key2: best peptides
my %Peptide_Copy_Number=();

open (my $IN,"<",$Motifs_Hits_File) || die "can't open IN '$Motifs_Hits_File' $!";
my $Current_motif_name="NA";
while (my $line=<$IN>)
{
	chomp ($line);
	if ($line=~/### Total (\d+) hits of PSSM '(\S+)'/) # new Motif
	{
		$Current_motif_name=$2;
		$Motifs_Name_To_Header{$Current_motif_name}=$line;
		push (@motifs,$Current_motif_name);
	}
	else # hit line
	{
		#NAHPLSTS        90.7068 850_Length_8_Repeats_90.7068208393131_Type_8    -63.8314        -4      ----NAHPLSTS----
#		YAPSSVGPSC      0.241242        33939_Length_10_Repeats_0.241241544785407_Type_10       -107.094        -2      --YAPSSVGPSC-------------
		my ($hit_peptide,$copy_number,$peptide_name,$score)=split(/\t/,$line);
		$hit_peptide=~s/^\s+|\s+$//g;
		if ($hit_peptide eq "YAPSSVGPSC")
		{
			my $wait=1;
		}
		if (exists $Motif_Best_Hits{$hit_peptide}{BEST_SCORE})
		{
			if ($score>$Motif_Best_Hits{$hit_peptide}{BEST_SCORE}) # update
			{
				# delete the hit from the previous motifs assigned to
				delete $Motivs_hits{$Motif_Best_Hits{$hit_peptide}{BEST_MOTIVE}}{$hit_peptide};
				# update
				$Motif_Best_Hits{$hit_peptide}{BEST_SCORE}=$score;
				$Motif_Best_Hits{$hit_peptide}{BEST_MOTIVE}=$Current_motif_name;
				$Motivs_hits{$Current_motif_name}{$hit_peptide}=$line;
			}
		}
		else # assign
		{
			$Motif_Best_Hits{$hit_peptide}{BEST_SCORE}=$score;
			$Motif_Best_Hits{$hit_peptide}{BEST_MOTIVE}=$Current_motif_name;
			$Motivs_hits{$Current_motif_name}{$hit_peptide}=$line;
		}
		if (!exists $Peptide_Copy_Number{$hit_peptide})
		{
			$Peptide_Copy_Number{$hit_peptide}=$copy_number;
		}
	}
}
close ($IN);

open (my $OUT,">",$Out_File) || die "Can't open OUT '$Out_File' $!";
foreach my $motif (@motifs)
{
	my $Num_Of_Unique_hits=0;
	if (exists $Motivs_hits{$motif})
	{
		$Num_Of_Unique_hits=scalar (keys $Motivs_hits{$motif});
	}
	$Motifs_Name_To_Header{$motif}=~s/Total \d+ hits/Total $Num_Of_Unique_hits hits/; # update the number of motifs hits
	print $OUT "$Motifs_Name_To_Header{$motif}\n";
	if (exists $Motivs_hits{$motif})
	{
		foreach my $hit (sort {$Peptide_Copy_Number{$b}<=>$Peptide_Copy_Number{$a}} keys $Motivs_hits{$motif})
		{
			print $OUT "$Motivs_hits{$motif}{$hit}\n";
		}
	}
}
close ($OUT);