#!/usr/bin/perl

open (OUTFILE, ">$ARGV[0]/gencode.qduplications.tsv")
	or die "Cannot open gencode.qduplications.tsv\n";

open (BENCHFILE, "$ARGV[0]/gencode.pc_annotations.out.tsv")
	or die "Cannot open gencode.pc_annotations.out.tsv\n";

open (PRINFILE, "$ARGV[0]/appris.principals.tsv")
	or die "Cannot open appris.principals.tsv\n";

%principal = ();
while (<PRINFILE>)
	{
	chomp;
	@id = split /\t/;
	$prin = $id[2];
	$principal{$id[1]} = $id[2];
	}


@trans = ();
%seq = ();
%annot = ();
$old_gene = "";
while (<BENCHFILE>)
	{
	chomp;
	@id = split /\t/;
	$gene = $id[0];
	$tran = $id[1];
	$annot = $id[4];
	if ($old_gene ne $gene)
		{
		for ($a=0;$a<=$#trans;$a++)
			{
			$tag = 0;
# print "$trans[$a]\t$principal{$gene}\n";
			if ($trans[$a] eq $principal{$old_gene} && ($annot{$trans[$a]} eq "NMD" || $annot{$trans[$a]} eq "NSD" || $annot{$trans[$a]} eq "PoP" || $annot{$trans[$a]} eq "RT"))
				{ print OUTFILE "$old_gene\t$trans[$a]\tPrincipal.$annot{$trans[$a]}\t$annot{$trans[$a]}\n"; $tag = 1}
			elsif ($trans[$a] eq $principal{$old_gene})
				{ print OUTFILE "$old_gene\t$trans[$a]\tPrincipal\t$annot{$trans[$a]}\n"; $tag = 1}
			elsif ($seq{$principal{$old_gene}} eq $seq{$trans[$a]})
				{ print OUTFILE "$old_gene\t$trans[$a]\tPrincipal Duplication\t$annot{$trans[$a]}\n"; $tag = 1}
			elsif ($seq{$principal{$old_gene}} =~ $seq{$trans[$a]} && ($annot{$trans[$a]} eq "SNF" || $annot{$trans[$a]} eq "ENF"))
				{ print OUTFILE "$old_gene\t$trans[$a]\tRedundant Principal\t$annot{$trans[$a]}\n"; $tag = 1}
			elsif ($annot{$trans[$a]} eq "NMD" || $annot{$trans[$a]} eq "NSD" || $annot{$trans[$a]} eq "PoP" || $annot{$trans[$a]} eq "RT")
				{ print OUTFILE "$old_gene\t$trans[$a]\tAlternative.$annot{$trans[$a]}\t$annot{$trans[$a]}\n"; $tag = 1}
			else
				{
				$count = $a+1;
				for ($p=$count;$p<=$#trans;$p++)
					{
					if ($seq{$trans[$a]} eq $seq{$trans[$p]})
						{
						print OUTFILE "$old_gene\t$trans[$a]\tAlternative Duplication\t$annot{$trans[$a]}\n";
						$tag = 1;
						last;
						}
					elsif ($seq{$trans[$p]} =~ $seq{$trans[$a]} && ($annot{$trans[$a]} eq "SNF" || $annot{$trans[$a]} eq "ENF"))
						{
						print OUTFILE "$old_gene\t$trans[$a]\tRedundant Alternative\t$annot{$trans[$a]}\t$trans[$p]\n";
						$tag = 1;
						last;
						}
					}
				}
			if (!$tag)
				{ print OUTFILE "$old_gene\t$trans[$a]\tAlternative\t$annot{$trans[$a]}\n"; }
			}
		@trans = ();
		}
	push @trans, $tran;
	$seq{$tran} = $id[2];
	$annot{$tran} = $annot;
	$old_gene = $gene;
	}


for ($a=0;$a<=$#trans;$a++)
	{
	$tag = 0;
# print "$trans[$a]\t$principal{$old_gene}\t$seq{$trans[$a]}\t$seq{$principal{$old_gene}}\t$annot{$trans[$a]}\n";
	if ($trans[$a] eq $principal{$old_gene} && ($annot{$trans[$a]} eq "NMD" || $annot{$trans[$a]} eq "NSD" || $annot{$trans[$a]} eq "PoP" || $annot{$trans[$a]} eq "RT"))
		{ print OUTFILE "$old_gene\t$trans[$a]\tPrincipal.$annot\t$annot{$trans[$a]}\n"; $tag = 1}
	elsif ($trans[$a] eq $principal{$old_gene})
		{ print OUTFILE "$old_gene\t$trans[$a]\tPrincipal\t$annot{$trans[$a]}\n"; $tag = 1}
	elsif ($seq{$principal{$old_gene}} eq $seq{$trans[$a]})
		{ print OUTFILE "$old_gene\t$trans[$a]\tPrincipal Duplication\t$annot{$trans[$a]}\n"; $tag = 1}
	elsif ($seq{$principal{$old_gene}} =~ $seq{$trans[$a]} && ($annot{$trans[$a]} eq "SNF" || $annot{$trans[$a]} eq "ENF"))
		{ print OUTFILE "$old_gene\t$trans[$a]\tRedundant Principal\t$annot{$trans[$a]}\n"; $tag = 1}
	elsif ($annot{$trans[$a]} eq "NMD" || $annot{$trans[$a]} eq "NSD" || $annot{$trans[$a]} eq "PoP" || $annot{$trans[$a]} eq "RT")
		{ print OUTFILE "$old_gene\t$trans[$a]\tAlternative.$annot\t$annot{$trans[$a]}\n"; $tag = 1}
	else
		{
		$count = $a+1;
		for ($p=$count;$p<=$#trans;$p++)
			{
			if ($seq{$trans[$a]} eq $seq{$trans[$p]})
				{
				print OUTFILE "$old_gene\t$trans[$a]\tAlternative Duplication\t$annot{$trans[$a]}\n";
				$tag = 1;
				last;
				}
			elsif ($seq{$trans[$p]} =~ $seq{$trans[$a]} && ($annot{$trans[$a]} eq "SNF" || $annot{$trans[$a]} eq "ENF"))
				{
				print OUTFILE "$old_gene\t$trans[$a]\tRedundant Alternative\t$annot{$trans[$a]}\t$trans[$p]\n";
				$tag = 1;
				last;
				}
			}
		}
	if (!$tag)
		{ print OUTFILE "$old_gene\t$trans[$a]\tAlternative\t$annot{$trans[$a]}\n"; }
	}
