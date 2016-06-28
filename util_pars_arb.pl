#!/usr/bin/perl -w 

use strict;

my %store=("gg_id"=>1, "ncbi_acc_w_ver"=>1, "ncbi_gi"=>1,"isolation_source"=>1,"authors"=>1,"title"=>1,"greengenes_tax_string"=>1,"aligned_seq"=>1);
my @Fields=qw {gg_id ncbi_acc_w_ver ncbi_gi isolation_source authors title greengenes_tax_string aligned_seq};
open (OUT, "> gg_13_5_arb_summary");
print OUT join("\t", @Fields), "\tseq_length\n";
for(my $i=0; $i<=126; $i++)
{
	print $i, "\n";
	my $inf="gg_13_5_arb_".$i.".txt";
	open (IN, "$inf");
	while(defined(my $text=<IN>))
	{
		chomp $text;
		if($text eq "BEGIN")
		{
			my %store1=();
			foreach my $field (keys %store)
			{
				$store1{$field}=$store{$field};
				
			}
	
			while($text ne "END")
			{
				$text=<IN>;
				chomp $text;
				if($text ne "END")
				{
					$text=~/(.*)=(.*)/;
					my $field=$1;
					my $value=$2;
					#print $text,"\n";
					if(defined($store1{$field}))
					{
						if($store1{$field} eq "1")
						{
							if($field eq "aligned_seq")
							{
								$value=~s/-//g;
							}
							if($field eq "greengenes_tax_string")
							{
								if($value=~/^k__Bacteria/){$value="Bacteria";}
								elsif($value=~/^k__Archaea/){$value="Archaea";}
								else{ print $i, "--",$store1{"gg_id"}, "--", $value, "\n"; }
							}
							$store1{$field}=$value;
						}
						else
						{
							 print $i, "--",$field, "--", $store1{$field}, "--",$store1{"gg_id"},"\n";
						}
					}
				}
			}
			if(($store1{"isolation_source"} ne "") && ($store1{"authors"} ne "") && ($store1{"title"} ne ""))
			{
				foreach my $field ( @Fields)
				{
					print OUT $store1{$field}, "\t";
				}
				print OUT length($store1{"aligned_seq"}), "\n";
			}
		}
		
	}
	close(IN);
}
close(OUT);
