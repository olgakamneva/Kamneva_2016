#!/usr/bin/perl -w 

use strict;
my $tab_in="gg_13_5_arb_summary";
open (INT, $tab_in);
open (OUTT, "> gg_13_5_arb_summary_event_OTU");
while(defined(my $text1=<INT>))
{
	chomp $text1;
	my @a=split(/\t/, $text1);
	print OUTT $a[0], "\t", $a[3], "\t", $a[4], "\t", $a[5], "\n";
	
}
close (INT);
close (OUTT);

