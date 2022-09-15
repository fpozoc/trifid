#!/usr/bin/perl

open (PRINFILE, "$ARGV[0]/spade_references.tsv")
	or die "Cannot open pfam_reference.tsv\n";

open (NOTFILE, "$ARGV[0]/pc_translations.tsv")
	or die "Cannot open benchfile pc_translations.tsv\n";

open (OUT2FILE, ">$ARGV[0]/ListOfSawnOffTranscripts.tsv")
	or die "Cannot outfile\n";

open (OUTFILE, ">$ARGV[1]")
	or die "Cannot outfile\n";



$old_gene = "";
@transcripts = ();
@starts = ();
@ends = ();
@seqs = ();
@alternatives = ();
$damage = 0;
$genedamage = 0;
$genelost = 0;
$lost = 0;
$total_events = 0;
$gene_events = 0;

%hgnc = ();
%principal = ();
%principal1 = ();
@principals = ();
while (<PRINFILE>)
	{
	chomp;
	$_ =~ s/\"//g;
	@ids = split /\t/;
	$gene = $ids[1];
#print "$gene\n";
	# next if ($gene =~ /^ENSGR/);
        if ($old_gene && ($gene ne $old_gene))
                { 
		$principal{$old_gene} = join ";", @principals;
		$principal1{$old_gene} = $principals[0];
		$hgnc{$old_gene} = $old_hgnc;
#print "Here's $principal1{$old_gene} and $old_gene\n";
		@principals = ();
		}
	push @principals, $ids[2];
	$old_gene = $gene;
	$old_hgnc = $ids[0];
	}

$principal{$old_gene} = join ";", @principals;
$hgnc{$old_gene} = $old_hgnc;

$old_gene = 0;
$old_gene_suff = 0;
%alternative = ();
@alternatives = ();
while (<NOTFILE>)
	{
#exit if ($count > 10);
	chomp;
	$_ =~ s/\"//g;
	@ids = split /\t/;
	$gene = $ids[0];
	@gene = split /\./, $gene;
	# next if ($gene =~ /^ENSGR/);
	$transcript = $ids[1];
	@tranid = split /\./, $transcript;
        if ($old_gene && ($gene[0] ne $old_gene))
                { 
		$alternative{$old_gene} = join ";", @alternatives;
		@alternatives = ();
		}
	@principals = split ";", $principal{$gene};
	$principal = 0;
	for ($a=0;$a<=$#principals;$a++)
		{
		if ($tranid[0] eq $principals[$a])
			{$principal = 1}
		}
	if (!$principal)
		{ push @alternatives, $tranid[0];}
	$old_gene = $gene[0];
	}

$alternative{$old_gene} = join ";", @alternatives;

seek NOTFILE, 0, 0;

while (<NOTFILE>)
	{
#exit if ($count > 10);
	chomp;
	$_ =~ s/\"//g;
	@ids = split /\t/;
	$gene = $ids[0];
	@gene = split /\./, $gene; 
#exit if ($gene[0] eq "ENSG00000000460");
#if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")
 #if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")         
#if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")
 #if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")         
#if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")
 #if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")         
#if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")
 #if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")         
#if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")
 #if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")         
#if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")
 #if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")         
#if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")
 #if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")         
#if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")
 #if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")         
#if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")
 #if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")         
#if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")
 #if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")         
#if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")
 #if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")         
#if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")
 #if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")         
#if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")
 #if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")         
#if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")
 #if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")         
#if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")
 #if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")         
#if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")
 #if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")         
#if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")
 #if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")         
#if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")
 #if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")         
#if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")
 #if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")         
#if ($gene[0] ne "ENSG00000130402" && $gene[0] ne "ENSG00000140416" && $gene[0] ne "ENSG00000169714" && $gene[0] ne "ENSG00000122367")
#if ($gene[0] ne "ENSG00000183091" && $gene[0] ne "ENSG00000125741")
 #if ($gene[0] ne "ENSG00000183091" && $gene[0] ne "ENSG00000125741")                                                                           
#if ($gene[0] ne "ENSG00000183091" && $gene[0] ne "ENSG00000125741")
 #if ($gene[0] ne "ENSG00000183091" && $gene[0] ne "ENSG00000125741")                                                                           
#if ($gene[0] ne "ENSG00000183091" && $gene[0] ne "ENSG00000125741")
 #if ($gene[0] ne "ENSG00000183091" && $gene[0] ne "ENSG00000125741")                                                                           
#if ($gene[0] ne "ENSG00000183091" && $gene[0] ne "ENSG00000125741")
 #if ($gene[0] ne "ENSG00000183091" && $gene[0] ne "ENSG00000125741")                                                                           
#if ($gene[0] ne "ENSG00000183091" && $gene[0] ne "ENSG00000125741")
 #if ($gene[0] ne "ENSG00000183091" && $gene[0] ne "ENSG00000125741")                                                                           
#if ($gene[0] ne "ENSG00000183091" && $gene[0] ne "ENSG00000125741")
 #if ($gene[0] ne "ENSG00000183091" && $gene[0] ne "ENSG00000125741")                                                                           
#if ($gene[0] ne "ENSG00000183091" && $gene[0] ne "ENSG00000125741")
 #if ($gene[0] ne "ENSG00000183091" && $gene[0] ne "ENSG00000125741")                                                                           
#if ($gene[0] ne "ENSG00000183091" && $gene[0] ne "ENSG00000125741")
 #if ($gene[0] ne "ENSG00000183091" && $gene[0] ne "ENSG00000125741")                                                                           
#if ($gene[0] ne "ENSG00000183091" && $gene[0] ne "ENSG00000125741")
 #if ($gene[0] ne "ENSG00000183091" && $gene[0] ne "ENSG00000125741")                                                                           
#if ($gene[0] ne "ENSG00000183091" && $gene[0] ne "ENSG00000125741")
 #if ($gene[0] ne "ENSG00000183091" && $gene[0] ne "ENSG00000125741")                                                                           
#if ($gene[0] ne "ENSG00000183091" && $gene[0] ne "ENSG00000125741")
 #if ($gene[0] ne "ENSG00000183091" && $gene[0] ne "ENSG00000125741")                                                                           
#if ($gene[0] ne "ENSG00000183091" && $gene[0] ne "ENSG00000125741")
     #   if ($gene[0] ne "ENSG00000013573" && $gene[0] ne "ENSG00000142168" && $gene[0] ne "ENSG00000178209" && $gene[0] ne "ENSG00000122367"    && $gene[0] ne "ENSG00000183091")                                                                                                              
     #   {next}                                                                                                                                 
	# next if ($gene =~ /^ENSGR/);
 #   next if ($gene !~ /ENSG00000077380/ && $gene !~ /ENSG00000004866/ && $gene !~ /ENSG00000128595/ && $gene !~ /ENSG00000204628/);            
	$transcript = $ids[1];
	@tranid = split /\./, $transcript;
	$seq = $ids[2];
	
	chop $seq;
        if ($old_gene && ($gene[0] ne $old_gene))
                {
					# print $principal1{$old_gene};
		if (!open (SPADFILE, "$ARGV[0]/spade3/$principal1{$old_gene}.pfam"))
			{ 
				# print "Cannot open spadefile $principal1{$old_gene} at $old_gene\n"; 
			}
# print "$old_gene and principal $principal1{$old_gene}\n";
		$gene_events = 0;
		$genedamage = 0;
		$genelost = 0;
		$pfam = 0;
		while (<SPADFILE>)
			{
			@stuff = split /\t/;
			$start = $stuff[2];
			$end = $stuff[3];
			$pfam = $stuff[1];
			push @starts, $start;
			push @ends, $end;
			push @pfamnames, $pfam;
#print "$start $end $pfam\n";
			}

		@alternatives = split ";", $alternative{$old_gene};
#print "This is $alternative{$old_gene}\n";
#print "And this is $principal1{$old_gene} and $old_gene\n";
		for ($a=0;$a<=$#alternatives;$a++)
			{
			@nevents = ();
			@events = ();
			@cevents = ();
#print "Thru\n";
			if (open (PAIRFILE, "$ARGV[0]/muscle/$old_gene.$principal1{$old_gene}.$alternatives[$a].ali"))
				{}
			else
				{print "$old_gene.$principal1{$old_gene}.$alternatives[$a].ali fail!\n"}
			$transcriptpair = join " ", $alternatives[$a], $principal{$old_gene};
#print "$transcriptpair\n";
			&analyse_align;
	
			my %found = ();
			my @events_u = grep { !$found{$_}++ } @events; 
			my %found = ();
			my @cevents_u = grep { !$found{$_}++ } @cevents; 
			my %found = ();
			my @nevents_u = grep { !$found{$_}++ } @nevents; 
#if ($old_gene eq "ENSG00000183091")
#		{print "$#nevents and $#nevents_u vs $#events and $#events_u vs $#cevents and $#cevents_u\n"}
#print "$#nevents and $#nevents_u vs $#events and $#events_u vs $#cevents and $#cevents_u\n";
			$lost = 0;
			$tranevents = 0;
			$damage = 0;
			for ($d=0;$d<=$#events_u;$d++)
				{
				@resevent = split "_", $events_u[$d];
				$delete = $resevent[1] - $resevent[0];
				$subst = $delete - $resevent[6];
				if ($resevent[2] && $delete < 2)
					{$reslost = $resevent[2]; $type = "Insert"}
				elsif ($subst < 5)
					{$reslost = $delete; $type = "Delete"}
				else
					{$reslost = $delete; $type = "Subst"}
				$homonogaps = $resevent[5] - $resevent[4];
				if ($resevent[7])
					{$type = "NAGNAG"; }
				elsif ($resevent[3] > 30 && $resevent[4] <= 10 && $reslost > 10)
					{$type = "Homology"; }
				elsif ($resevent[3] >= 30 && $resevent[6] <= 5 && $reslost > 24)
					{$type = "Homology"; }
				elsif ($resevent[3] >= 30 && $resevent[6] <= 2 && $reslost > 11)
					{$type = "Homology"; }
				elsif ($resevent[5] >= 40 && $resevent[4] <= 10 && $reslost > 10)
					{$type = "Homology"; }
				elsif ($resevent[5] >= 40 && $resevent[6] <= 4 && $reslost > 24)
					{$type = "Homology"; }
				elsif ($resevent[5] >= 40 && $resevent[6] <= 1 && $reslost > 11)
					{$type = "Homology"; }
				elsif ($resevent[5] >= 40 && !$resevent[6] && $reslost > 9)
					{$type = "Homology"; }
				elsif ($homonogaps > 20 && $reslost > 9)
					{$type = "Homology"; }
				elsif ($homonogaps >= 50 && $reslost < 11 && $reslost >= 8)
					{$type = "Homology"; }
				elsif ($homonogaps >= 55 && $reslost < 8 && $reslost >= 5)
					{$type = "Homology"; }
#if 5 or less gaps and iden > 30
#if 5 or less gaps and iden > 30
#if 15% gaps and iden > 30
#if homology > 33 and 0 gaps
#if homology > 40 and 15% gaps or 5 gaps
				for ($p=0;$p<=$#starts;$p++)
					{
					$loss = 0;
					$pfamlost = 0;
					$pfamdisturb = 0;
					$damages = 0;
#print "$alternatives[$a]\t$resevent[0], $starts[$p], $resevent[1], $ends[$p]\n";
					if ($resevent[0] <= $starts[$p] && $resevent[1] >= $ends[$p])
						{
						$loss=1; 
						$pfamlost = $ends[$p] - $starts[$p];
						$pfamdisturb = $pfamlost;
						}
					elsif ($resevent[0] > $starts[$p] && $resevent[0] < $ends[$p] && $resevent[1] > $starts[$p] && $resevent[1] < $ends[$p] && $type eq "Insert")
						{
						if ($reslost > 12)
							{
							$damages=1;
							$pfamlost = $reslost;
							$pfamdisturb = $ends[$p] - $starts[$p];
							}
						elsif ($reslost > 4)
							{
							$damages=1;
							$pfamlost = $reslost;
							$pfamdisturb = ($ends[$p] - $starts[$p])/2;
							}
						else
							{
							$damages=0;
							$pfamlost = $reslost;
							$pfamdisturb = ($ends[$p] - $starts[$p])/10;
							}
						}
					elsif ($resevent[0] > $starts[$p] && $resevent[0] < $ends[$p] && $resevent[1] > $starts[$p] && $resevent[1] < $ends[$p])
						{
						$diff1 = $resevent[0] - $starts[$p];
						$diff2 = $ends[$p] - $resevent[1];
						if ($diff1 < 4 && $diff2 < 4)
							{$loss=1;}
						elsif ($delete > 4)
							{$damages=1;}
						elsif ($type eq "Insert" && $reslost > 4)
							{$damages=1;}
						if ($reslost > $delete)
							{ $pfamlost = $reslost; }
						else
							{ $pfamlost = $delete; }
						$pfamdisturb = $pfamlost;
						}
					elsif ($resevent[0] > $starts[$p] && $resevent[0] < $ends[$p])
						{
						$diff1 = $resevent[0] - $starts[$p];
						$diff2 = $ends[$p] - $resevent[0];
						if ($diff1 < 4)
							{$loss=1;}
						elsif ($diff2 > 4)
							{$damages=1;}
						$pfamlost = $ends[$p] - $resevent[0];
						$pfamdisturb = $pfamlost;
						}
					elsif ($resevent[1] > $starts[$p] && $resevent[1] < $ends[$p])
						{
						$diff1 = $resevent[1] - $starts[$p];
						$diff2 = $ends[$p] - $resevent[1];
						if ($diff2 < 4)
							{$loss=1;}
						elsif ($diff1 > 4)
							{$damages=1;}
						$pfamlost = $resevent[1] - $starts[$p];
						$pfamdisturb = $pfamlost;
						}
					if ($loss)
						{$result = "Lost"}
					elsif ($damages)
						{$result = "Damaged"}
					else
						{$result = "Intact"}
					if ($type eq "Homology")
						{$damages = 0; $lost = 0}
					if ($resevent[0] <= $starts[$p] && $resevent[1] >= $ends[$p] && !$damages)
						{ $loss = 1; }
					$pc_lost = $pfamdisturb/($ends[$p]-$starts[$p])*100;
					print OUTFILE "$old_gene\t$transcriptpair\t$p\t$pfamnames[$p]\t$starts[$p]\t$ends[$p]\t$type\t$reslost\t$pfamlost\t$result\t$resevent[3]\t$resevent[4]\t$pc_lost\t$hgnc{$old_gene}\n";
					}
				if (!$pfam)
					{ print OUTFILE "$old_gene\t$transcriptpair\tNULL\tNULL\tNULL\tNULL\t$type\t$reslost\tNULL\tNULL\t$resevent[3]\t$resevent[4]\t\t$hgnc{$old_gene}\n"; }
					
				if ($damages)
					{$damage++;$genedamage++}
#print "$damages and $damage\n";
				if ($loss)
					{$lost++;$genelost++}
				$total_events++;
				$gene_events++;
				}
			$tranevents = $d;
#print "$alternatives[$a] at indel at $lost and $damage and $tranevents\n";
			for ($d=0;$d<=$#cevents_u;$d++)
				{
				@resevent = split "_", $cevents_u[$d];
				$reslost = $resevent[1] - $resevent[0];
				if ($resevent[4] == 100)
					{$type = "C-term Delete"}
				else
					{$type = "C-term Swap"}
				$homonogaps = $resevent[5] - $resevent[4];
if ($alternatives[$a] eq "ENST00000358278")
	{print "$resevent[3]\t$resevent[4]\t$resevent[5]\t$resevent[6]\t$homonogaps\t$reslost\n";}
				if ($resevent[3] >= 30 && $resevent[4] <= 10 && $reslost > 10)
					{$type = "Homology"; }
				elsif ($resevent[3] >= 30 && $resevent[6] <= 5 && $reslost > 24)
					{$type = "Homology"; }
				elsif ($resevent[3] >= 30 && $resevent[6] <= 2 && $reslost > 11)
					{$type = "Homology"; }
				elsif ($resevent[5] >= 40 && $resevent[4] <= 10 && $reslost > 10)
					{$type = "Homology"; }
				elsif ($resevent[5] >= 40 && $resevent[6] <= 4 && $reslost > 24)
					{$type = "Homology"; }
				elsif ($resevent[5] >= 40 && $resevent[6] <= 1 && $reslost > 11)
					{$type = "Homology"; }
				elsif ($resevent[5] >= 40 && !$resevent[6] && $reslost > 9)
					{$type = "Homology"; }
				elsif ($homonogaps > 33 && $reslost > 9)
					{$type = "Homology"; }
				elsif ($homonogaps >= 50 && $reslost < 11 && $reslost >= 8)
					{$type = "Homology"; }
				elsif ($homonogaps >= 55 && $reslost < 8 && $reslost >= 5)
					{$type = "Homology"; }
#if 5 or less gaps and iden > 30
#if 15% gaps and iden > 30
#if homology > 33 and 0 gaps
#if homology > 40 and 15% gaps or 5 gaps
				for ($p=0;$p<=$#starts;$p++)
					{
#print "@resevent\n";
					$loss = 0;
					$damages = 0;
					$pfamlost = 0;
					$pfamdisturb = 0;
#print "$alternatives[$a]\t$resevent[0], $starts[$p], $resevent[1], $ends[$p]\n";
					if ($resevent[0] > $starts[$p] && $resevent[0] < $ends[$p])
						{
						$diff1 = $resevent[0] - $starts[$p];
						$diff2 = $ends[$p] - $resevent[0];
						if ($diff1 < 5)
							{$loss=1; }
						elsif ($diff2 > 4)
							{$damages=1; }
						$pfamlost = $ends[$p] - $resevent[0];
						}
					elsif ($resevent[0] <= $starts[$p])
						{ $loss = 1; $pfamlost = $ends[$p] - $starts[$p]}
					$pfamdisturb = $pfamlost;
					if ($loss)
						{$result = "Lost";}
					elsif ($damages)
						{$result = "Damaged"; }
					else
						{$result = "Intact"; }
					if ($type eq "Homology")
						{$damages = 0; $lost = 0}
					if ($resevent[2] eq "Two Proteins")
						{print OUTFILE "$old_gene\t$transcriptpair\t$p\t$pfamnames[$p]\t$starts[$p]\t$ends[$p]\tTwo Proteins\t$reslost\t$pfamlost\t$result\t$resevent[3]\t$resevent[4]\t\t$hgnc{$old_gene}\n"; }
					else
						{
						$pc_lost = $pfamdisturb/($ends[$p]-$starts[$p])*100;
						print OUTFILE "$old_gene\t$transcriptpair\t$p\t$pfamnames[$p]\t$starts[$p]\t$ends[$p]\t$type\t$reslost\t$pfamlost\t$result\t$resevent[3]\t$resevent[4]\t$pc_lost\t$hgnc{$old_gene}\n"; 
						}
					}
				if (!$pfam && $resevent[2] ne "Two Proteins")
					{
					print OUTFILE "$old_gene\t$transcriptpair\tNULL\tNULL\tNULL\tNULL\t$type\t$reslost\tNULL\tNULL\t$resevent[3]\t$resevent[4]\t\t$hgnc{$old_gene}\n";
					}

				elsif (!$pfam && $resevent[2] eq "Two Proteins")
					{print OUTFILE "$old_gene\t$transcriptpair\tNULL\tNULL\tNULL\tNULL\tTwo Proteins\t$reslost\tNULL\tNULL\t$resevent[3]\t$resevent[4]\t\t$hgnc{$old_gene}\n";}
				if ($damages)
					{$damage++;$genedamage++}
				if ($loss)
					{$lost++;$genelost++}
				$total_events++;
				$gene_events++;
				}
			$tranevents = $tranevents+$d;
#print "$alternatives[$a] at cterm at $lost and $damage and $tranevents\n";
			for ($d=0;$d<=$#nevents_u;$d++)
				{
				@resevent = split "_", $nevents_u[$d];
				$reslost = $resevent[0];
				if (!$resevent[2])
					{$type = "N-term Delete"}
				else
					{$type = "N-term Swap"}
				$homonogaps = $resevent[4] - $resevent[3];
				if ($resevent[2] > 30 && $resevent[3] <= 10 && $reslost > 10)
					{$type = "Homology"; }
				elsif ($resevent[2] >= 30 && $resevent[5] <= 5 && $reslost > 24)
					{$type = "Homology"; }
				elsif ($resevent[2] >= 30 && $resevent[5] <= 2 && $reslost > 11)
					{$type = "Homology"; }
				elsif ($resevent[4] >= 40 && $resevent[3] <= 10 && $reslost > 10)
					{$type = "Homology"; }
				elsif ($resevent[4] >= 40 && $resevent[5] <= 4 && $reslost > 24)
					{$type = "Homology"; }
				elsif ($resevent[4] >= 40 && $resevent[5] <= 1 && $reslost > 11)
					{$type = "Homology"; }
				elsif ($resevent[4] >= 40 && !$resevent[5] && $reslost > 9)
					{$type = "Homology"; }
				elsif ($homonogaps > 33 && $reslost > 9)
					{$type = "Homology"; }
				elsif ($homonogaps >= 50 && $reslost < 11 && $reslost >= 8)
					{$type = "Homology"; }
				elsif ($homonogaps >= 55 && $reslost < 8 && $reslost >= 5)
					{$type = "Homology"; }
				for ($p=0;$p<=$#starts;$p++)
					{
					$loss = 0;
					$damages = 0;
					$pfamlost = 0;
					$pfamdisturb = 0;
					if ($resevent[0] > $starts[$p] && $resevent[0] < $ends[$p])
						{
						$diff1 = $resevent[0] - $starts[$p];
						$diff2 = $ends[$p] - $resevent[0];
						if ($diff2 < 5)
							{$loss=1; }
						elsif ($diff1 > 4)
							{$damages=1; }
						$pfamlost = $resevent[0] - $starts[$p];
						}
					elsif ($resevent[0] >= $ends[$p])
						{ $loss = 1; $pfamlost = $ends[$p] - $starts[$p]}
					$pfamdisturb = $pfamlost;
					if ($loss)
						{$result = "Lost"; }
					elsif ($damages)
						{$result = "Damaged"; }
					else
						{$result = "Intact"}
					if ($type eq "Homology")
						{$damages = 0; $lost = 0}
					$pc_lost = $pfamdisturb/($ends[$p]-$starts[$p])*100;
					print OUTFILE "$old_gene\t$transcriptpair\t$p\t$pfamnames[$p]\t$starts[$p]\t$ends[$p]\t$type\t$resevent[0]\t$pfamlost\t$result\t$resevent[2]\t$resevent[3]\t$pc_lost\t$hgnc{$old_gene}\n";
					}
				if (!$pfam)
					{print OUTFILE "$old_gene\t$transcriptpair\tNULL\tNULL\tNULL\tNULL\t$type\t$resevent[0]\tNULL\tNULL\t$resevent[2]\t$resevent[3]\t\t$hgnc{$old_gene}\n";}
				if ($damages)
					{$damage++;$genedamage++}
				if ($loss)
					{$lost++;$genelost++}
				$total_events++;
				$gene_events++;
				}
			$tranevents = $tranevents+$d;
##print "$alternatives[$a] at nterm at $lost and $damage and $tranevents\n";
			}

		@transcripts = ();
		@starts = ();
		@pfamnames = ();
		@ends = ();
		@seqs = ();
		@alternatives = ();
		print OUTFILE "$old_gene\tLost = $genelost\tDamaged = $genedamage\tEvents = $gene_events\n";
		}
	push @transcripts, $tranid[0];
	push @seqs, $seq;
	$old_gene = $gene[0];
	}
	


#print "$old_gene and principal $principal1{$old_gene}\n";
$gene_events = 0;
$genedamage = 0;
$genelost = 0;
$pfam = 0;
if (!open (SPADFILE, "$ARGV[0]/spade3/$principal1{$old_gene}.pfam"))
	{
		# print "Cannot open spadefile $principal1{$old_gene} at $old_gene\n"; 
		}
else
	{
	while (<SPADFILE>)
		{
		@stuff = split /\t/;
		$start = $stuff[2];
		$end = $stuff[3];
		$pfam = $stuff[1];
		push @starts, $start;
		push @ends, $end;
		push @pfamnames, $pfam;
#print "$old_gene $principal1{$old_gene} $start $end $pfam\n";
		}

	@alternatives = split ";", $alternative{$old_gene};
#print "$alternative{$old_gene}\n";
	for ($a=0;$a<=$#alternatives;$a++)
		{
		@events = ();
		@nevents = ();
		@cevents = ();
#print "Thru\n";
		if (open (PAIRFILE, "$ARGV[0]/muscle/$old_gene.$principal1{$old_gene}.$alternatives[$a].ali"))
			{}
		else
			{print "$old_gene.$principal1{$old_gene}.$alternatives[$a].ali fail!"}
		$transcriptpair = join " ", $alternatives[$a], $principal{$old_gene};
#print "$transcriptpair\n";
		&analyse_align;

	
		my %found = ();
		my @events_u = grep { !$found{$_}++ } @events; 
		my %found = ();
		my @cevents_u = grep { !$found{$_}++ } @cevents; 
		my %found = ();
		my @nevents_u = grep { !$found{$_}++ } @nevents; 
			
#print "$#nevents and $#nevents_u vs $#events and $#events_u vs $#cevents and $#cevents_u\n";
		$damage = 0;
		$lost = 0;
		$tranevents = 0;
		$damage = 0;
		for ($d=0;$d<=$#events_u;$d++)
			{
			$tranevents = 0;
			@resevent = split "_", $events_u[$d];
#print "Resevent: $resevent[0],$starts[0] and $resevent[1],$ends[0]\n";
			$delete = $resevent[1] - $resevent[0];
			$subst = $delete - $resevent[6];
			if ($resevent[2] && $delete < 2)
				{$reslost = $resevent[2]; $type = "Insert"}
			elsif ($subst < 5)
				{$reslost = $delete; $type = "Delete"}
			$homonogaps = $resevent[5] - $resevent[4];
			if ($resevent[7])
				{$type = "NAGNAG"; }
			elsif ($resevent[3] > 30 && $resevent[4] <= 10 && $reslost > 10)
				{$type = "Homology"; }
			elsif ($resevent[3] >= 30 && $resevent[6] <= 5 && $reslost > 24)
				{$type = "Homology"; }
			elsif ($resevent[3] >= 30 && $resevent[6] <= 2 && $reslost > 11)
				{$type = "Homology"; }
			elsif ($resevent[5] >= 40 && $resevent[4] <= 10 && $reslost > 10)
				{$type = "Homology"; }
			elsif ($resevent[5] >= 40 && $resevent[6] <= 4 && $reslost > 24)
				{$type = "Homology"; }
			elsif ($resevent[5] >= 40 && $resevent[6] <= 1 && $reslost > 11)
				{$type = "Homology"; }
			elsif ($resevent[5] >= 40 && !$resevent[6] && $reslost > 9)
				{$type = "Homology"; }
			elsif ($homonogaps > 33 && $reslost > 9)
				{$type = "Homology"; }
			elsif ($homonogaps >= 50 && $reslost < 11 && $reslost >= 8)
				{$type = "Homology"; }
			elsif ($homonogaps >= 55 && $reslost < 8 && $reslost >= 5)
				{$type = "Homology"; }
			for ($p=0;$p<=$#starts;$p++)
				{
				$loss = 0;
				$damages = 0;
				$pfamlost = 0;
				$pfamdisturb = 0;
#print "$alternatives[$a]\t$resevent[0], $starts[$p], $resevent[1], $ends[$p]\n";
				if ($resevent[0] <= $starts[$p] && $resevent[1] >= $ends[$p])
					{
					$loss=1; 
					$pfamlost = $ends[$p] - $starts[$p];
					$pfamdisturb = $pfamlost;
					}
				elsif ($resevent[0] >= $starts[$p] && $resevent[0] <= $ends[$p] && $resevent[1] >= $starts[$p] && $resevent[1] <= $ends[$p] && $type eq "Insert")
					{
					if ($reslost > 12)
						{
						$damages=1;
						$pfamdisturb = $ends[$p] - $starts[$p];
						$pfamlost = $reslost;
						}
					elsif ($reslost > 4)
						{
						$damages=1;
						$pfamdisturb = ($ends[$p] - $starts[$p])/2;
						$pfamlost = $reslost;
						}
					else
						{
						$damages=0;
						$pfamdisturb = ($ends[$p] - $starts[$p])/10;
						$pfamlost = $reslost;
						}
					}
				elsif ($resevent[0] >= $starts[$p] && $resevent[0] <= $ends[$p] && $resevent[1] >= $starts[$p] && $resevent[1] <= $ends[$p])
					{
					$diff1 = $resevent[0] - $starts[$p];
					$diff2 = $ends[$p] - $resevent[1];
					if ($diff1 < 4 && $diff2 < 4)
						{$loss=1;}
					elsif ($delete > 4)
						{$damages=1;}
					$pfamlost = $delete;
					$pfamdisturb = $pfamlost;
					}
				elsif ($resevent[0] >= $starts[$p] && $resevent[0] <= $ends[$p])
					{
					$diff1 = $resevent[0] - $starts[$p];
					$diff2 = $ends[$p] - $resevent[0];
					if ($diff1 < 4)
						{$loss=1;}
					elsif ($diff2 > 4)
						{$damages=1;}
					$pfamlost = $ends[$p] - $resevent[0];
					$pfamdisturb = $pfamlost;
					}
				elsif ($resevent[1] >= $starts[$p] && $resevent[1] <= $ends[$p])
					{
					$diff1 = $resevent[1] - $starts[$p];
					$diff2 = $ends[$p] - $resevent[1];
					if ($diff2 < 4)
						{$loss=1;}
					elsif ($diff1 > 4)
						{$damages=1;}
					$pfamlost = $resevent[1] - $starts[$p];
					$pfamdisturb = $pfamlost;
					}
#print "Indel: $diff1 and $diff2, $delete, $resevent[2] \n";
				if ($resevent[0] <= $starts[$p] && $resevent[1] >= $ends[$p] && !$damages)
					{ $loss = 1; }
				if ($loss)
					{$result = "Lost"}
				elsif ($damages)
					{$result = "Damaged"}
				else
					{$result = "Intact"}
#print "For Indels: $resevent[2] and $delete and $resevent[3] and $damages\n";
				if ($type eq "Homology")
					{$damages = 0; $lost = 0}

				$pc_lost = $pfamdisturb/($ends[$p]-$starts[$p])*100;
				print OUTFILE "$old_gene\t$transcriptpair\t$p\t$pfamnames[$p]\t$starts[$p]\t$ends[$p]\t$type\t$reslost\t$pfamlost\t$result\t$resevent[3]\t$resevent[4]\t$pc_lost\t$hgnc{$old_gene}\n";
				}
			if (!$pfam)
				{print OUTFILE "$old_gene\t$transcriptpair\tNULL\tNULL\tNULL\tNULL\t$type\t$reslost\tNULL\tNULL\t$resevent[3]\t$resevent[4]\t\t$hgnc{$old_gene}\n";}
			if ($damages)
				{$damage++;$genedamage++}
#print "$damages and $damage\n";
			if ($loss)
				{$lost++;$genelost++}
			$total_events++;
			$gene_events++;
			}
		$tranevents = $d;
#print "$alternatives[$a] at indel at $lost and $damage and $tranevents\n";
		for ($d=0;$d<=$#cevents_u;$d++)
			{
			@resevent = split "_", $cevents_u[$d];
			$reslost = $resevent[1] - $resevent[0];
			if ($resevent[4] == 100)
				{$type = "C-term Delete"}
			else
				{$type = "C-term Swap"}
			$homonogaps = $resevent[5] - $resevent[4];
			if ($resevent[3] > 30 && $resevent[4] <= 10 && $reslost > 10)
				{$type = "Homology"; }
			elsif ($resevent[3] >= 30 && $resevent[6] <= 5 && $reslost > 24)
				{$type = "Homology"; }
			elsif ($resevent[3] >= 30 && $resevent[6] <= 2 && $reslost > 11)
				{$type = "Homology"; }
			elsif ($resevent[5] >= 40 && $resevent[4] <= 10 && $reslost > 10)
				{$type = "Homology"; }
			elsif ($resevent[5] >= 40 && $resevent[6] <= 4 && $reslost > 24)
				{$type = "Homology"; }
			elsif ($resevent[5] >= 40 && $resevent[6] <= 1 && $reslost > 11)
				{$type = "Homology"; }
			elsif ($resevent[5] >= 40 && !$resevent[6] && $reslost > 9)
				{$type = "Homology"; }
			elsif ($homonogaps > 33 && $reslost > 9)
				{$type = "Homology"; }
			elsif ($homonogaps >= 50 && $reslost < 11 && $reslost >= 8)
				{$type = "Homology"; }
			elsif ($homonogaps >= 55 && $reslost < 8 && $reslost >= 5)
				{$type = "Homology"; }
			for ($p=0;$p<=$#starts;$p++)
				{
				$loss = 0;
				$damages = 0;
				$pfamlost = 0;
				$pfamdisturb = 0;
				if ($resevent[0] > $starts[$p] && $resevent[0] < $ends[$p])
					{
					$diff1 = $resevent[0] - $starts[$p];
					$diff2 = $ends[$p] - $resevent[0];
					if ($diff1 < 5)
						{$loss=1; }
					elsif ($diff2 > 4)
						{$damages=1; }
					$pfamlost = $ends[$p] - $resevent[0];
					}
				elsif ($resevent[0] <= $starts[$p])
					{ $loss = 1; $pfamlost = $ends[$p] - $starts[$p]}
				$pfamdisturb = $pfamlost;
				if ($loss)
					{$result = "Lost"; }
				elsif ($damages)
					{$result = "Damaged"; }
				else
					{$result = "Intact"; }
				if ($type eq "Homology")
					{$damages = 0; $lost = 0}
				if ($resevent[2] eq "Two Proteins")
					{print OUTFILE "$old_gene\t$transcriptpair\t$p\t$pfamnames[$p]\t$starts[$p]\t$ends[$p]\tTwo Proteins\t$reslost\t$pfamlost\t$result\t$resevent[3]\t$resevent[4]\t\t$hgnc{$old_gene}\n"; }
				else
					{
					$pc_lost = $pfamdisturb/($ends[$p]-$starts[$p])*100;
					print OUTFILE "$old_gene\t$transcriptpair\t$p\t$pfamnames[$p]\t$starts[$p]\t$ends[$p]\t$type\t$reslost\t$pfamlost\t$result\t$resevent[3]\t$resevent[4]\t$pc_lost\t$hgnc{$old_gene}\n"; 
					}
			if (!$pfam && $resevent[2] ne "Two Proteins")
				{ print OUTFILE "$old_gene\t$transcriptpair\tNULL\tNULL\tNULL\tNULL\t$type\t$reslost\tNULL\tNULL\t$resevent[3]\t$resevent[4]\t\t$hgnc{$old_gene}\n"; }
			elsif (!$pfam && $resevent[2] eq "Two Proteins")
				{print OUTFILE "$old_gene\t$transcriptpair\tNULL\tNULL\tNULL\tNULL\tTwo Proteins\t$reslost\tNULL\tNULL\t$resevent[3]\t$resevent[4]\t\t$hgnc{$old_gene}\n";}
			if ($damages)
				{$damage++;$genedamage++}
			if ($loss)
				{$lost++;$genelost++}
			$total_events++;
			$gene_events++;
			}
		$tranevents = $tranevents+$d;
#print "$alternatives[$a] at cterm at $lost and $damage and $tranevents\n";
		for ($d=0;$d<=$#nevents_u;$d++)
			{
			@resevent = split "_", $nevents_u[$d];
			$reslost = $resevent[0];
			if (!$resevent[2])
				{$type = "N-term Delete"}
			else
				{$type = "N-term Swap"}
			$homonogaps = $resevent[4] - $resevent[3];
			if ($resevent[2] > 30 && $resevent[3] <= 10 && $reslost > 10)
				{$type = "Homology"; }
			elsif ($resevent[2] >= 30 && $resevent[5] <= 5 && $reslost > 24)
				{$type = "Homology"; }
			elsif ($resevent[2] >= 30 && $resevent[5] <= 2 && $reslost > 11)
				{$type = "Homology"; }
			elsif ($resevent[4] >= 40 && $resevent[3] <= 10 && $reslost > 10)
				{$type = "Homology"; }
			elsif ($resevent[4] >= 40 && $resevent[5] <= 4 && $reslost > 24)
				{$type = "Homology"; }
			elsif ($resevent[4] >= 40 && $resevent[5] <= 1 && $reslost > 11)
				{$type = "Homology"; }
			elsif ($resevent[4] >= 40 && !$resevent[5] && $reslost > 9)
				{$type = "Homology"; }
			elsif ($homonogaps > 33 && $reslost > 9)
				{$type = "Homology"; }
			elsif ($homonogaps >= 50 && $reslost < 11 && $reslost >= 8)
				{$type = "Homology"; }
			elsif ($homonogaps >= 55 && $reslost < 8 && $reslost >= 5)
				{$type = "Homology"; }
			for ($p=0;$p<=$#starts;$p++)
				{
				$loss = 0;
				$damages = 0;
				$pfamlost = 0;
				$pfamdisturb = 0;
				if ($resevent[0] > $starts[$p] && $resevent[0] < $ends[$p])
					{
					$diff1 = $resevent[0] - $starts[$p];
					$diff2 = $ends[$p] - $resevent[0];
					if ($diff2 < 5)
						{$loss=1; }
					elsif ($diff1 > 4)
						{$damages=1; }
					$pfamlost = $resevent[0] - $starts[$p];
					}
				elsif ($resevent[0] >= $ends[$p])
					{ $loss = 1; $pfamlost = $ends[$p] - $starts[$p]}
				$pfamdisturb = $pfamlost;
				if ($loss)
					{$result = "Lost";}
				elsif ($damages)
					{$result = "Damaged"; }
				else
					{$result = "Intact"}
				if ($type eq "Homology")
					{$damages = 0; $lost = 0}
				$pc_lost = $pfamdisturb/($ends[$p]-$starts[$p])*100;
				print OUTFILE "$old_gene\t$transcriptpair\t$p\t$pfamnames[$p]\t$starts[$p]\t$ends[$p]\t$type\t$resevent[0]\t$pfamlost\t$result\t$resevent[2]\t$resevent[3]\t$pc_lost\t$hgnc{$old_gene}\n";
				}
			if (!$pfam)
				{print OUTFILE "$old_gene\t$transcriptpair\tNULL\tNULL\tNULL\tNULL\t$type\t$resevent[0]\tNULL\tNULL\t$resevent[2]\t$resevent[3]\t\t$hgnc{$old_gene}\n";}
			if ($damages)
				{$damage++;$genedamage++}
			if ($loss)
				{$lost++;$genelost++}
			$total_events++;
			$gene_events++;
			}
		$tranevents = $tranevents+$d;
#print "$alternatives[$a] at nterm at $lost and $damage and $tranevents\n";
		}
	}

print OUTFILE "$old_gene\tLost = $genelost\tDamaged = $genedamage\tEvents = $gene_events\n";
print OUTFILE "Lost = $genelost\nDamaged = $genedamage\nAll = $total_events\n";




sub analyse_align
        {
#print "Thru sub at $old_gene\n";
        @seq1 = ();
        @seq2 = ();
        $string = <PAIRFILE>;
        $string = <PAIRFILE>;
        $string = <PAIRFILE>;
        $count = 0;
        while (<PAIRFILE>)
                {
                chomp;
                if (!$count)
                        { @data1 = split " "; push @seq1, $data1[2]}
                elsif ($count == 1)
                        { @data2 = split " "; push @seq2, $data2[2]}
                $count++;
                if ($count > 3)
                        {$count = 0}
                }
        $seq1 = join ("", @seq1);
        $seq1 =~ s/X/-/g;
        $seq2 = join ("", @seq2);
        $seq2 =~ s/X/-/g;
#print "$seq1, $seq2\n";
        if ($seq1 eq $seq2)
                { next; }
	$prinseq = $seq1;
	$prinseq =~ s/-//g;
	$prinlen = length($prinseq);
        @seq1 = split "", $seq1;
	@seq2 = split "", $seq2;
#print "$seq1\n";

	$pc_iden=0;
	$pc_gaps=0;
	$homo = 0;
	$iden = 0;
	$cterm = 0;
	$gaps = 0;
	$subtract = 0;
	$seqno = $prinlen;
	$not_same_c = 0;
	$cterm_iden=0;
	$cterm_gaps=0;
        for ($m=$#seq1;$m>0;$m--)
                {
#print "Cterm $m, $cterm, $seqno, $iden, $cterm_iden,$seq1[$m], $seq2[$m]\n";
                if ($seq1[$m] ne $seq2[$m] && $seq1[$m] ne "-" && $seq2[$m] ne "-")
                        {
			$cterm = $seqno; 
			$iden = 0;
			if ($seq1[$m] =~ /[FILMV]/ && $seq2[$m] =~ /[FILMV]/)
				{$homo++}
			if ($seq1[$m] =~ /[ILM]/ && $seq2[$m] =~ /[ILM]/)
				{$homo++}
			if ($seq1[$m] =~ /[YF]/ && $seq2[$m] =~ /[YF]/)
				{$homo++}
			if ($seq1[$m] =~ /[HQ]/ && $seq2[$m] =~ /[HQ]/)
				{$homo++}
			if ($seq1[$m] =~ /[ND]/ && $seq2[$m] =~ /[ND]/)
				{$homo++}
			if ($seq1[$m] =~ /[ED]/ && $seq2[$m] =~ /[ED]/)
				{$homo++}
			if ($seq1[$m] =~ /[EQ]/ && $seq2[$m] =~ /[EQ]/)
				{$homo++}
			if ($seq1[$m] =~ /[HRK]/ && $seq2[$m] =~ /[HRK]/)
				{$homo++}
			$not_same_c = 1
			}
                elsif ($seq1[$m] eq $seq2[$m] && $seq1[$m] ne "-")
                        {
			$cterm_iden++;
			$iden++
			}
		if ($seq1[$m] eq "-" || $seq2[$m] eq "-")
			{
			$cterm = $seqno; 
			$gaps++;
			$iden = 0;
			if ($seq2[$m] eq "-")
				{ $cterm_gaps++; }
			}
		if ($iden > 3 && !$cterm)
			{ last }
		elsif ($iden > 3 && !$not_same_c)
			{
			if ($gaps and !$not_same_c)
				{print OUT2FILE "$transcriptpair\tC\n"}
#NOTE GET RID OF THIS COMMENT IF YOU WANT TO COUNT REAL EVENTS			last
                	}
		if ($iden > 6)
			{ last }
		if ($seq1[$m] ne "-")
			{$seqno--}
                }
#print "This is $m in Two Proteins\n";
	if ($m)
		{ $pause_c = $m;}
	else
		{ $pause_c = $m}
	$c_diff = $#seq1-$m-$iden+1;
	if ($c_diff)
		{ $pc_iden = ($cterm_iden- $iden)*100/$c_diff; }
	if ($c_diff)
		{ $pc_homo = ($cterm_iden-$iden+($homo/2))*100/$c_diff; }
	if ($gaps >= $c_diff)
		{ $pc_gaps = 100; }
	else
		{ 
		if ($c_diff)
			{ $pc_gaps = $gaps*100/$c_diff; }
		}
	if (!$m && $cterm < 3)
		{
		$type = "Two Proteins";
		$cevent = join "_", $cterm, $prinlen, $type, $pc_iden, $pc_gaps, $pc_homo, $gaps;
#print "$alternatives[$a]\tCterm\t$cevent\n";
		push @cevents, $cevent;
		}
	else
		{
		if ($cterm ne "0")
			{
#print "This is Ciden: $cterm_iden vs $iden\n";
#			$cterm--;
			$cevent = join "_", $cterm, $prinlen, $c_diff, $pc_iden, $pc_gaps, $pc_homo, $gaps; 
#print "$alternatives[$a]\tCterm\t$cevent\n";
			push @cevents, $cevent;
			}
		}
#print "$alternatives[$a]\t$cterm\n";

	$homo = 0;
	$iden = 0;
	$gaps = 0;
	$seqno = 1;
	$nterm = 0;
	$not_same_n = 0;
	$pc_iden=0;
	$pc_gaps=0;
	$nterm_iden=0;
	$nterm_gaps=0;
        for ($m=0;$m<$pause_c;$m++)
                {
                if ($seq1[$m] ne $seq2[$m] && $seq1[$m] ne "-" && $seq2[$m] ne "-")
                        {
			$nterm = $seqno; 
			$iden = 0;
			if ($seq1[$m] =~ /[FILMV]/ && $seq2[$m] =~ /[FILMV]/)
				{$homo++}
			if ($seq1[$m] =~ /[ILM]/ && $seq2[$m] =~ /[ILM]/)
				{$homo++}
			if ($seq1[$m] =~ /[ND]/ && $seq2[$m] =~ /[ND]/)
				{$homo++}
			if ($seq1[$m] =~ /[ED]/ && $seq2[$m] =~ /[ED]/)
				{$homo++}
			if ($seq1[$m] =~ /[YF]/ && $seq2[$m] =~ /[YF]/)
				{$homo++}
			if ($seq1[$m] =~ /[HQ]/ && $seq2[$m] =~ /[HQ]/)
				{$homo++}
			if ($seq1[$m] =~ /[EQ]/ && $seq2[$m] =~ /[EQ]/)
				{$homo++}
			if ($seq1[$m] =~ /[HRK]/ && $seq2[$m] =~ /[HRK]/)
				{$homo++}
			$not_same_n = 1
			}
                elsif ($seq1[$m] eq $seq2[$m] && $seq1[$m] ne "-")
                        {
			$nterm_iden++; 
			$iden++
			}
		if ($seq1[$m] eq "-" || $seq2[$m] eq "-")
			{
			$nterm = $seqno; 
			$gaps++;
			$iden = 0;
			if ($seq2[$m] eq "-")
				{ $nterm_gaps++; }
                	}
		if ($iden > 3 && !$nterm)
			{ last }
		if ($iden > 3 && !$not_same_n)
			{
			if ($gaps and !$not_same_n)
				{print OUT2FILE "$transcriptpair\tN\n"}
#NOTE GET RID OF THIS COMMENT IF YOU WANT TO COUNT REAL EVENTS			last
                	}
		if ($iden > 6)
			{ last }

#print "Nterm $nterm, $seqno, $m, $iden, $gaps, $nterm_gaps, $seq1[$m] - $seq2[$m]\n";
		if ($seq1[$m] ne "-")
			{$seqno++}
#print "stuck d!\n";
                }
	if ($nterm)
		{ 
		$n_diff = $m-$iden+1;
		$pc_iden = ($nterm_iden-$iden)*100/$n_diff;
		$pc_homo = ($nterm_iden-$iden+($homo/2))*100/$n_diff; 
		if ($gaps >= $nterm)
			{ $pc_gaps = 100; }
		else
			{ $pc_gaps = $gaps*100/$n_diff; }
		$nevent = join "_", $nterm, $n_diff, $pc_iden, $pc_gaps, $pc_homo, $gaps; 
#print "$nevent\n";
#print "$alternatives[$a]\tNterm\t$nevent\n";
		push @nevents, $nevent; 
		}
##print "$alternatives[$a]\t$nterm\n";

#print "$seqno B4\n";
#	$pause_n = $m;
	if (!$nterm)
		{$pause_n = 0}
	else
		{$pause_n = $m - 3}
#	$seqno = $nterm;
	$homo = 0;
	$iden = 0;
	$internvent = 0;
	$insertion = 0;
	$seqno = $seqno - 4;
#	$seqno = $seqno;
#print "$seqno after\n";
	$event_start=0;
	$event_end = 0;
	$event_iden=0;
	$event_gaps=0;
	$pc_iden=0;
	$pc_gaps=0;
#print "Thru $pause_n, $pause_c, $seqno\n";
        for ($m=$pause_n;$m<=$pause_c;$m++)
                {
#print "Swap $swap_start, $swap_end, $seqno, $m, $iden, $swap_iden, $swap_gaps, $seq1[$m], $seq2[$m]\n";
#print "Indel $indel_start, $indel_end, $seqno, $m, $iden, $indel_iden, $indel_gaps, $seq1[$m], $seq2[$m]\n";
                if ($seq1[$m] ne $seq2[$m] && $seq1[$m] ne "-" && $seq2[$m] ne "-")
                        {
			if (!$event_start)
				{$event_start = $seqno}
			if ($seq1[$m] =~ /[FILMV]/ && $seq2[$m] =~ /[FILMV]/)
				{$homo++}
			if ($seq1[$m] =~ /[ILM]/ && $seq2[$m] =~ /[ILM]/)
				{$homo++}
			if ($seq1[$m] =~ /[YF]/ && $seq2[$m] =~ /[YF]/)
				{$homo++}
			if ($seq1[$m] =~ /[HQ]/ && $seq2[$m] =~ /[HQ]/)
				{$homo++}
			if ($seq1[$m] =~ /[ND]/ && $seq2[$m] =~ /[ND]/)
				{$homo++}
			if ($seq1[$m] =~ /[ED]/ && $seq2[$m] =~ /[ED]/)
				{$homo++}
			if ($seq1[$m] =~ /[EQ]/ && $seq2[$m] =~ /[EQ]/)
				{$homo++}
			if ($seq1[$m] =~ /[HRK]/ && $seq2[$m] =~ /[HRK]/)
				{$homo++}
			$event_end = $seqno; 
			$iden=0;
			}
                elsif ($seq1[$m] eq $seq2[$m] && $seq1[$m] ne "-")
                        {
			if ($event_start)
				{ $event_iden++; }
			$iden++
			}
		elsif ($seq1[$m] eq "-" || $seq2[$m] eq "-")
			{
			if (!$event_start)
				{ $event_start = $seqno; }
			$event_end = $seqno; 
			$event_gaps++; 
			if ($seq1[$m] eq "-" && $seq2[$m] ne "-")
				{$insertion++}
			$iden=0;
#print "$m, $seqno, $iden, $insertion, $indel_start, $swap_start\n";
			}
		if ($seq1[$m] ne "-")
			{$seqno++; }
		if ($iden > 10 && $event_start)
			{
			$event_end++;
			$event_len = $event_end-$event_start;
			if ($event_len < 5)
				{$pc_iden = 0; $pc_homo = 0}
			else
				{$pc_iden = ($event_iden-$iden)*100/$event_len; $pc_homo = ($event_iden-$iden+($homo/2))*100/$event_len;}
			if ($event_gaps >= $event_len)
				{ $pc_gaps = 100; }
			else
				{ $pc_gaps = $event_gaps*100/$event_len; }
			$NN_tag = $event_len - $event_gaps;
#print "This one: $indel_gaps, $insertion, $indel_len\n";
			if ($event_len < 5 && $event_gaps && $event_gaps < 5)
				{$NN_tag = 1}
			elsif (($event_gaps == $insertion && $event_gaps < 5 && $event_len < 2) || ($event_gaps == $insertion && $event_gaps < 6 && $event_len == 2)) 
				{$NN_tag = 1}
			else
				{$NN_tag = 0}			$internevent = join "_", $event_start, $event_end, $insertion, $pc_iden, $pc_gaps, $pc_homo, $event_gaps, $NN_tag; 
#print "$alternatives[$a]\tSWAP\t$internvent\n";
			if ($internevent ne "0")
				{push @events, $internevent;} 

			$event_start=0;
			$event_end=0;
			$event_iden=0;
			$event_gaps=0;
			$insertion=0;
			$iden = 0;
			$homo = 0;
			}
		}
	if ($event_start)
		{
		$event_end++;
		$event_len = $event_end-$event_start;
		if ($event_len < 5)
			{$pc_iden = 0; $pc_homo = 0}
		else
			{$pc_iden = ($event_iden-$iden)*100/$event_len; $pc_homo = ($event_iden-$iden+($homo/2))*100/$event_len;}
		if ($event_gaps >= $event_len)
			{ $pc_gaps = 100; }
		else
			{ $pc_gaps = $event_gaps*100/$event_len; }
		if ($event_len < 5 && $event_gaps && $event_gaps < 5)
			{$NN_tag = 1}
		elsif (($event_gaps == $insertion && $event_gaps < 5 && $event_len < 2) || ($event_gaps == $insertion && $event_gaps < 6 && $event_len == 2)) 
			{$NN_tag = 1}
		else
			{$NN_tag = 0}
		$internevent = join "_", $event_start, $event_end, $insertion, $pc_iden, $pc_gaps, $pc_homo, $event_gaps, $NN_tag; 
#print "$alternatives[$a]\tSWAP\t$internevent\n";
		if ($internevent ne "0")
			{push @events, $internevent;} 

		}
	}
	}

