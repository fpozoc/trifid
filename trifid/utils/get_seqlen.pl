#!/usr/bin/perl

open (OUTFILE, ">$ARGV[0]/gencode.pc_annotations.out.tsv")
	or die "Cannot open gencode.pc_annotations.out.tsv\n";

open (ANNOTFILE, "$ARGV[0]/gencode.pc_annotations.tsv")
	or die "Cannot open gencode.pc_annotations.tsv\n";

open (BENCHFILE, "$ARGV[0]/appris.pc_sequences.tsv")
	or die "Cannot open appris.pc_sequences.tsv\n";

%annot = ();
while (<ANNOTFILE>)
	{
	chomp;
	$annot = "";
	@id = split /\t/;
    if ($id[4] eq "nonsense_mediated_decay")
		{$annot = NMD}
	elsif ($id[4] eq "non_stop_decay")
		{$annot = NSD}
	elsif ($id[4] eq "polymorphic_pseudogene")
		{$annot = PoP}
	elsif ($id[5] eq "readthrough")
		{$annot = RT}
	elsif ($id[6] eq "End NF")
		{$annot = ENF}
	elsif ($id[6] eq "Start NF")
		{$annot = SNF}
	$annot{$id[2]} = $annot;
	}

while (<BENCHFILE>)
	{
	chomp;
	@id = split /\t/;
	@gene = $id[0];
	@data = $id[1];
	$id[2] =~ s/Z//g;
	$id[2] =~ s/X//g;
	$seqlen = length($id[2]);
	print OUTFILE "$gene[0]\t$data[0]\t$id[2]\t$seqlen\t$annot{$data[0]}\n";
	}
