#!/usr/bin/perl -s
#
#quick and dirty converter for paired end read data in interleaved format to
#two separate files containing matching reads. optionally, also renames reads
#to create matching read names in the two output files.
#
#usage: interleaved_to_1_2.pl [-renamed] interleaved.fastq

use warnings;
use Getopt::Std;

open IN, "<$ARGV[0]" or die "cannot open input file $ARGV[0]: $!";
$ARGV[0]=~s/(.*)\.fastq/$1/; 
open OUT1, ">$ARGV[0].1.fastq" or die "cannot open first output file: $!";
open OUT2, ">$ARGV[0].2.fastq" or die "cannot open second output file: $!";

if (!$renamed) {
	while (<IN>){
		print OUT1 $_;
		my $next=<IN>; print OUT1 $next;
		$next=<IN>; print OUT1 $next;
		$next=<IN>; print OUT1 $next;
		$next=<IN>; print OUT2 $next;
		$next=<IN>; print OUT2 $next;
		$next=<IN>; print OUT2 $next;
		$next=<IN>; print OUT2 $next;
	}
}
else {
	while (<IN>){
		$_=~s/(.*)\/.*/$1/; 
		my $name1 = $_;
                my $next2=<IN>;
		my $next3=<IN>;
		my $next4=<IN>;
                my $name2=<IN>;
		$name5=~s/(.*)\/.*/$1/;
		if ($name1 ne $name5) { next; }
                my $next6=<IN>;
                my $next7=<IN>;
                my $next8=<IN>;
		print OUT1 $name1, $next2, $next3, $next4;
		print OUT2 $name5, $next6, $next7, $next8;
        }
}

