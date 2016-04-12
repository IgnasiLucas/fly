#!/usr/bin/perl
#
# Cx001.plx
#
# This is meant to in silico cut a genome with one or two restriction sites.
# I will use Parallel::ForkManager to cut more than one scaffold or chromosome
# at a time. The output I want is only the distrubution of length fragments,
# but an option could be offered to produce the sequences. To deal with Ns in
# the reference genome, I envision two options. First, I can assume that there
# is no cut site within a tract of Ns. Otherwise, I can eliminate the fragments
# with Ns from the statistics.
#
# Output: fragment length, comment (OK|Ns|First|Last), number of masked bases
# in first 100 bases of the fragment, number of masked bases in last 100 bases
# of the fragment.

use lib '/home/ilucas/.lib/';
use warnings;
use strict;
use Getopt::Std;
use Parallel::ForkManager;
use Bio::SeqIO;

my %opts = ('m' => 2, 'f' => "fasta");
getopts('g:f:r:m:s', \%opts);
unless ((exists $opts{'g'}) && (exists $opts{'r'})) {
	print "Usage: Cx001.plx -g <genome.fasta> -r <restriction_enzyme_list> [-s] [-m]\n";
	print "   -g <genome.file>  \tSequence file of genome to be cut.\n";
	print "   -f <genome.format>\tFormat of the sequence file [fasta].\n";
	print "   -r <rest.enz.list>\tFile with list of restriction enzymes and combinations\n";
	print "                     \tthereof to be used. Text file with two columns: restriction\n";
	print "                     \tenzyme identifier, and sequence that it cuts. If more than\n";
	print "                     \tone enzyme is combined, a identify the combination and express\n";
	print "                     \tthe sequences cut as a perl pattern. The sequence is taken as a\n";
	print "                     \tstring representing a pattern.\n";
	print "   -s                \tWhether to output the sequences of the fragments or not [false].\n";
	print "   -m                \tMaximum number of parallel processes to run [2].\n";
	exit;
}

open(RES, $opts{'r'}) || die "I cannot open $opts{'r'}.\n";
my %restriction;
my $N = 0;
while (<RES>) {
	chomp;
	$N++;
	my @line = split /\t/, $_;
	if (scalar @line != 2) {
		die "Line $N of the restriction enzyme list file has more or less than 2 fields.\n";
	}
	$restriction{$line[0]} = $line[1];
	if (! -d $line[0]) {
		mkdir $line[0] or die "Error in mkdir $line[0]: $!\n";
	}
}
close RES;
my $manager = new Parallel::ForkManager( $opts{'m'} );
my $ChrNum = 0;
if ($N <= $opts{'m'}) { # if processing <= m restriction enzymes, they can be forked.
	for my $i (keys %restriction) {
		$manager->start and next;
		my $Genome = Bio::SeqIO->new(-file => $opts{'g'}, -format => $opts{'f'});
		my %Fragments;
		open(OUT, ">$i/lengths") or die "I cannot open $i/lengths.\n";
		if ($opts{'s'}) {
			open(SEQ, ">$i/fa") or die "I cannot open $i/seq.\n";
		}
		while (my $scaffold = $Genome->next_seq) {
			my $ScaffId = $scaffold->id;
			$Fragments{$ScaffId} = [ split /$restriction{$i}/i, $scaffold->seq ];
			my $fragId = 1;
			for my $f (@{ $Fragments{$ScaffId} }) {
				my $comment = "OK";
				my $masked1 = 0;
				my $end1 = substr($f, 0, 100);
				while ($end1 =~ /[acgtn]/g) { $masked1++ }
				my $masked2 = 0;
				my $end2 = substr($f, -100, 100);
				while ($end2 =~ /[acgtn]/g) { $masked2++ }
				if ($f =~ /N/i) {
					$comment = "Ns";
				}
				if ($fragId == 1) {
					$comment = $comment . "First";
				}
				if ($fragId == scalar @{ $Fragments{$ScaffId} }) {
					$comment = $comment . "Last";
				}
				print OUT length($f), "\t", $comment, "\t", $masked1, "\t", $masked2, "\t", $ScaffId, "\n";
				if ($opts{'s'}) { print SEQ ">$ScaffId.$fragId\n", $f, "\n" }
				$fragId++;
			}
		}
		close OUT;
		if ($opts{'s'}) { close SEQ }
		$manager->finish;
	}
} else { # Otherwise, I should fork the scaffolds or chromosomes.
	my $Genome = Bio::SeqIO->new(-file => $opts{'g'}, -format => $opts{'f'});
	while (my $scaffold = $Genome->next_seq) {
		$manager->start and next;
		my $ScaffId = $scaffold->id;
		for my $i (keys %restriction) {
			open(OUT, ">$i/$ScaffId.length") or print "I cannot open $i/$ScaffId.lengths.\n";
			if ($opts{'s'}) {
				open(SEQ, ">$i/$ScaffId.fa") or print "I cannot open $i/ScaffId.fa.\n";
			}
			my @Fragments = split /$restriction{$i}/i, $scaffold->seq;
			my $fragId = 1;
			for my $f (@Fragments) {
				my $comment = "OK";
				my $masked1 = 0;
                                my $end1 = substr($f, 0, 100);
                                while ($end1 =~ /[acgtn]/g) { $masked1++ }
                                my $masked2 = 0;
                                my $end2 = substr($f, -100, 100);
                                while ($end2 =~ /[acgtn]/g) { $masked2++ }
				if ($f =~ /N/i) {
					$comment = "Ns";
				}
				if ($fragId == 1) {
					$comment = $comment . "First";
				}
				if ($fragId == scalar @Fragments) {
					$comment = $comment . "Last";
				}
				print OUT length($f), "\t", $comment, "\t", $masked1, "\t", $masked2, "\t", $ScaffId, "\n";
				if ($opts{'s'}) { print SEQ ">$ScaffId.$fragId\n", $f, "\n" }
				$fragId++;
			}
			close OUT;
			if ($opts{'s'}) { close SEQ }
		}
		$manager->finish;
	}
}

$manager->wait_all_children;
