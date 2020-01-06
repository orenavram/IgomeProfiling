use strict;

if (@ARGV<2) {die "USAGE: perl $0 <Hits_File> <Out_Dir> <Motifs_Of_Interest_File?>\n";}
my $Hits_File=shift;
my $Out_Dir=shift;
my $Motifs_Of_Interest_File=shift;


if ($Out_Dir!~/\/$/){$Out_Dir=$Out_Dir."/";}
if (!-d $Out_Dir) {mkdir($Out_Dir);}

my %MotifsToExtract=();
if (defined $Motifs_Of_Interest_File)
{
	open (my $TO_EXTRACT,"<",$Motifs_Of_Interest_File) || die "Can't open TO_EXTRACT '$Motifs_Of_Interest_File' $!";
	while (my $line=<$TO_EXTRACT>)
	{
		chomp ($line);
		if (($line=~/,/) and ($Motifs_Of_Interest_File=~/.csv$/)) # To support RF best features format
		{
			my @line=split(/,/,$line);
			my $MotifName=$line[1];
			if ($MotifName=~/^"(\S+)"$/){$MotifName=$1;}
			$MotifsToExtract{$MotifName}=1;
		}
		else # regular list
		{
			$MotifsToExtract{$line}=1;
		}
	}
	close ($TO_EXTRACT);
}

open (my $HITS_FILE,"<",$Hits_File) || die "Can't open HITS_FILE '$Hits_File' $!";
my $line=<$HITS_FILE>;
while (defined $line)
{
	chomp ($line);
	if ($line=~/### Total (\d+) hits of PSSM '(\S+)'/) # new Motif
	{
		my $Number_of_hits=$1;
		if ($Number_of_hits>0)
		{
			my $Current_motif_name=$2;
			if ((!defined $Motifs_Of_Interest_File) or (exists $MotifsToExtract{$Current_motif_name}))
			{
				print "[INFO] Extract hits of motif $Current_motif_name\n";
				my $Current_motif_hits_file=$Out_Dir.$Current_motif_name.".hits.fas";
				$line=<$HITS_FILE>;
				open (my $FASTA,">",$Current_motif_hits_file) || die "Can't open FASTA '$Current_motif_hits_file' $!";
				while ((defined $line)&&($line!~/^###/)) # go over hits line: #NAHPLSTS        90.7068 850_Length_8_Repeats_90.7068208393131_Type_8    -63.8314        -4      ----NAHPLSTS----
				{
					my ($hit_peptide,$copy_number,$peptide_name,$score)=split(/\t/,$line);
					$hit_peptide=~s/^\s+|\s+$//g;
					print $FASTA ">$peptide_name\n$hit_peptide\n";
					$line=<$HITS_FILE>
				}
				close ($FASTA);
			}
			else
			{
				print "[INFO] Skipping motif: $Current_motif_name\n";
				$line=<$HITS_FILE>;
				while ((defined $line)&&($line!~/^###/)) # go over hits line but don't extract them
				{
					$line=<$HITS_FILE>;
				}
			}
		}
		else
		{
			print "[INFO] Skip motif: $Current_motif_name with 0 hits\n";
			$line=<$HITS_FILE>;
		}
	}
	else 
	{
		$line=<$HITS_FILE>;
	}
}
close ($HITS_FILE);