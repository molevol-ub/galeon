#!/home/vadim/miniconda3/envs/Galeon/bin/perl
use strict;
use warnings;

# Script to retrieve the list of scaffolds and its length from a fasta file

#usage: perl get_scaffold_length.pl genome_assembly.fasta

my ($name, $line, $nameout, $contig, $frame, $name2);
my (%nrfa, %inicio, %fin, %fastacutted, %fastanames);
my $skip = 0;
my $ids = "";

my $out = "";


if ($ARGV[0] =~ /(\S+).fasta/){
	$out = $1;
} else {
	$out = $ARGV[0];
}


my $file= "$ARGV[0]";
my %fasta;
open (Fasta , "<", $file);
while (<Fasta>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	if ($line =~ /^>(\S+)/) {
		$name = $1;
	} else {
		$fasta{$name} .= uc($line); # Change all sequences to upper case nucleotides
	}
}
close Fasta;

#open (Res, ">", "$ARGV[0]\_scaflonger100Kb.fasta");
open (Rest, ">", "$out\_scaffold_length.txt");
#open (Resren, ">", "$out\_final.fasta");

foreach my $key (sort keys %fasta){
	my $length = length ($fasta{$key});
	print Rest "$key\t$length\t$key\n";
	if ($length >= 100000){
#		print Res ">$key\n$fasta{$key}\n";
	}
}
#close Res;
close Rest;

# system ("sort -n -r -k 2 $out\_scaffold_length.txt > $out\_scaffold_length_sorted.txt");
system ("sort -n -r -k 2 $out\_scaffold_length.txt > ChrSizes.txt"); # update 24 Abril 2024
system ("rm $out\_scaffold_length.txt");

my $n = 1;

=skip
open (File , "<", "$out\_length_sorted.txt");
while (<File>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	my @subl = split (/\t/, $line);
	print Resren ">Scaffold$n\n$fasta{$subl[0]}\n";
	$n++;
}
close Fasta;
=cut



