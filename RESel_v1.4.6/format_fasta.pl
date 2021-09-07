#! /usr/bin/perl
use strict;
use warnings;

open OUT, ">$ARGV[1]" or die$!;
open INFILE, "$ARGV[0]" or die$!;

while(<INFILE>) {
	chomp($_);
	
	if ($_ =~ /^>/) {
			print OUT "$_\n";
	}
	
	else {
		print OUT "$_";
	}		
}	

close OUT;
close INFILE;
