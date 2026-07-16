#!/usr/bin/env/perl

# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# perl renameRMDLconsensi.pl <input> <prefix> <output>	#
# ===================================================== #
# Replaces sequence names from RepeatModeler standard	#
# output (consensi.fa and related files) with easy-to-	#
# use, short names + species identifier as prefix.	#
# ===================================================== #
# Alexander Suh				    12 Nov 2015 #
# Modification by Ricard Fontsere 07 Jul 2026#
#       Added compatibility to new RM2 outputs, where 
#    LTR pipeline consensi was giving duplicated names  #
# ===================================================== #
# Example:						#
# "perl renameRMDLconsensi.pl consensi.fa.classified \	#
# lepSin lepSin_rm1.0.lib"				#
# -----------------------------------------------------	#
# Old FASTA header: ">rnd-2_family-111#LINE/CR1-Zenon \	#
# ( Recon Family Size = 56, Final Multiple Alignment \	#
# Size = 53 )"						#
# -----------------------------------------------------	#
# New FASTA header: ">lepSin2-111#LINE/CR1-Zenon"	#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #

use strict;
use warnings;

my $file = $ARGV[0];
my $prefix = $ARGV[1];
my $out = $ARGV[2];

open (OUT, ">$out");

open (IN, $file);
while (my $line = <IN>) {
    my $typenum = "";
    if ($line =~ m/^>([^_]+)_/) {   # capture 'rnd-1' / 'ltr-2' from the raw header
	$typenum = $1;
	$typenum =~ s/-//g;          # -> 'rnd1' / 'ltr2'
    }
    $line =~ s/rnd-|ltr-/$prefix/g;
    $line =~ s/_family//g;
    if ($line =~ m/^>/) {
        chomp $line;
	my @parts = split(/ /, $line);
	my $id = $parts[0];
	if ($id =~ /#/) { $id =~ s/#/_$typenum#/; }   # insert before classification
	else            { $id .= "_$typenum"; }        # defensive: header with no '#'
	print OUT $id."\n";
    } else {
        printf OUT $line;
    }
}
close (IN);
close (OUT);