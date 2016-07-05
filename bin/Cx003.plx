#!/usr/bin/perl
#
# Here I implement an algorithm to write lists of 12 8-nucleotide
# words that I believe to be at a minimum Hamming distance of 4, and
# with balanced composition.

use warnings;
use strict;
use POSIX qw/floor/;

my $number = 12;
my $length = 8;
my @first;
my @alphabet = ("A", "C", "G", "T");
my %NextBase = ("A" => "C",
		"C" => "G",
		"G" => "T",
		"T" => "A",
		""  => "A" );
for (my $i = 0; $i < $number; $i++) {
	$first[$i] = $alphabet[$i % scalar @alphabet];
}

sub next_four {
	my @previous = split "", $_[0];
	my %index;
	for (my $i = 0; $i < $number; $i++) {
		push( @{ $index{$previous[$i]} }, $i);
	}
	my @options = ();
	for my $o (@alphabet) {
		my @current = ();
		$current[0] = $o;
		my $last = $o;
		for my $prev ($previous[0], $previous[1], $previous[2], $previous[3]) {
			for my $i (@{ $index{$prev} }) {
				unless (exists $current[$i]) {
					$current[$i] = $NextBase{$last};
					$last = $current[$i];
				}
			}
		}
		if (scalar @ current == $number) {
			push(@options, join("", @current));
		} else {
			print join("", @previous), "\n", join("-", @current), "\n";
			die "Error\n";
		}
	}
	return @options;
}

my %F = ();
my $stop = 0;
sub find_all {
	unless ((exists $F{$_[0]}) || ($stop == $length)) {
		push(@{ $F{$_[0]} }, next_four($_[0]));
	}
	for my $option (@{ $F{$_[0]} }) {
		unless ((exists $F{$option}) || ($stop == $length)) {
			find_all($option);
		}
	}
	$stop++;
}

my $first = join("", @first);
#print "# Building the tree of options.\n\n";
find_all($first);


for my $a (0, 1, 2, 3) {
	for my $b (0, 1, 2, 3) {
		for my $c (0, 1, 2, 3) {
			for my $d (0, 1, 2, 3) {
				for my $e (0, 1, 2, 3) {
					for my $f (0, 1, 2, 3) {
						for my $g (0, 1, 2, 3) {
							my @comb = (	                     $first,
									                  $F{$first}[$a],
									               $F{$F{$first}[$a]}[$b],
									            $F{$F{$F{$first}[$a]}[$b]}[$c],
									         $F{$F{$F{$F{$first}[$a]}[$b]}[$c]}[$d],
									      $F{$F{$F{$F{$F{$first}[$a]}[$b]}[$c]}[$d]}[$e],
									   $F{$F{$F{$F{$F{$F{$first}[$a]}[$b]}[$c]}[$d]}[$e]}[$f],
									$F{$F{$F{$F{$F{$F{$F{$first}[$a]}[$b]}[$c]}[$d]}[$e]}[$f]}[$g] );
							my @trans = ();
							for (my $i = 0; $i < $number; $i++) {
								$trans[$i] = substr($comb[0], $i, 1) . substr($comb[1], $i, 1) . substr($comb[2], $i, 1) . substr($comb[3], $i, 1) .
									     substr($comb[4], $i, 1) . substr($comb[5], $i, 1) . substr($comb[6], $i, 1) . substr($comb[7], $i, 1) ;
								print $trans[$i], "\t";
							}
							print "\n"
						}
					}
				}
			}
		}
	}
}
