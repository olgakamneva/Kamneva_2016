#!/usr/bin/perl -w 

# This splits species.mappings.v10.txt file from STRING database into 1 file per species
# and puts them into directory families in the current working directory, 
# if directory does not exist it will be created.
# To run: perl util_get_species_mappings.pl species.mappings.v10.txt



($inf1)=@ARGV;
print "Parsing species.mappings.v10.txt\n";

open (IN1, "$inf1");
$head=<IN1>;
%Maps=();
while(defined($text=<IN1>))
{
	chomp $text;
	my @a=split(/\s+/, $text);
	if(defined($Maps{$a[0]}))
	{
		push(@{$Maps{$a[0]}}, $text);
	}
	else
	{
		@{$Maps{$a[0]}}=();
		push(@{$Maps{$a[0]}}, $text);
	}
}
close(IN1);

print "Printing results\n";

unless(-e "families")
{
	mkdir "families";
}
foreach $sp (keys %Maps)
{
	$outf1="families/families_". $sp;
	open (OUT1, ">$outf1") || die "$!";
	print OUT1 $head;
	print OUT1 join("\n",@{$Maps{$sp}});
	close(OUT1);
		
}