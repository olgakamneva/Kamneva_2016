#!/usr/bin/perl -w 

use strict;
(my $seeds)=@ARGV;

my $inf="99_otu_map.txt";
open (IN, "$inf");
my %OTUs_99=();
my %OTUs_99_rev=();
while(defined(my $text=<IN>))
{
	chomp $text;
	my @a=split(/\t/, $text);
	my $otu=shift(@a);
	foreach my $o (@a)
	{
		${$OTUs_99{$otu}}{$o}=1;
		$OTUs_99_rev{$o}=$otu;
	}
}
close(IN);

print 26251, "\t", join("\t", keys %{$OTUs_99{26251}}),"\n";
print 533329, "\t", $OTUs_99_rev{533329},"\n";


$inf="97_otu_map.txt";
open (IN, "$inf");
my %OTUs_97=();
while(defined(my $text=<IN>))
{
	chomp $text;
	my @a=split(/\t/, $text);
	my $otu=shift(@a);
	foreach my $o (@a)
	{
		my $o99=$OTUs_99_rev{$o};
		foreach my $o99s (keys %{$OTUs_99{$o99}})
		{
			${$OTUs_97{$otu}}{$o99s}=1;
		}
	}
}
close(IN);



my $outf="otus/97_99_otu_map.txt";
open (OUT, "> $outf");

foreach my $o (keys %OTUs_97)
{
	foreach my $o99 (keys %{$OTUs_97{$o}})
	{
		print OUT $o, "\t", $o99, "\n";
	}
}
close(OUT);






