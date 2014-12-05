#!/uaopt/perl/5.14.2/bin/perl

# A quick converter of aligned reads in sam format to xml, for use with PIIM
# Assumes the input bam file has been sorted; otherwise may take extra time to sort while processing
# Based largely on examples provided by Lincoln Stein at http://search.cpan.org/~lds/Bio-SamTools/lib/Bio/DB/Sam.pm

use warnings;
use strict;
use Bio::DB::Sam;
use Bio::SeqIO;

# load bam file specified at the command line

my $sam = Bio::DB::Sam->new(-bam  =>"$ARGV[0]",
                             -fasta=>"$ARGV[1]");
# get length of reference genome
my $ref = Bio::SeqIO->new(-file => "$ARGV[1]", -format => "Fasta");
my $refseq = $ref->next_seq();
my $reflen = $refseq->length();
print STDERR $reflen;
							 
# print header information to outfile
print STDOUT "<?xml version=\"1.0\"?>
<!DOCTYPE Contig
SYSTEM \"piim_r.dtd\">
<Contig>", "\n";
							 
# load all pairs into one array? only those mapping to reference; still probably super memory-intensive
my @pairs = $sam->get_features_by_location(-type   => 'read_pair',
                                            -start  => 1,
                                            -end    => $reflen);

# introduce an array to keep track of coverage levels across the genome
my @coverage;

# for each pair, print required information in xml format (read start, read seq, read quality)
PAIRLOOP: for my $pair (@pairs) {
	my $length = $pair->length;   # insert length
	my ($first_mate,$second_mate) = $pair->get_SeqFeatures;
	if $first_mate->proper_pair {
		my $f_start = $first_mate->start;
		my $s_start = $second_mate->start;
		my $f_end = $first_mate->end;
		my $s_end = $second_mate->end;
		# noting coverage levels and skipping pair if over 50 at any bp
	        foreach my $i ($f_start..$f_end) {
			$coverage[$i]++;
			if ($coverage[$i]>$ARGV[2]) {
				next PAIRLOOP;
			}
		}
		foreach my $i ($s_start..$s_end) {
			$coverage[$i]++;
			if ($coverage[$i]>$ARGV[2]) {
				next PAIRLOOP;
			}
		}
		print STDERR $length, "\n";
		my $f_dna = $first_mate->query->dna;
		my $s_dna = $second_mate->query->dna;
		my @f_scores = $first_mate->qscore;
		my @s_scores = $second_mate->qscore;
		# converting qual scores
		for (my $i=0; $i < scalar(@f_scores); $i++) {
			$f_scores[$i]-=31;
		}
		for (my $j=0; $j < scalar(@s_scores); $j++) {
			$s_scores[$j]-=31;
		}
		print STDOUT " <Chr>
  <Read start=\"$f_start\">
   <Bases>$f_dna</Bases>
   <Quals>@f_scores</Quals>
  </Read>
  <Read start=\"$s_start\">
   <Bases>$s_dna</Bases>
   <Quals>@s_scores</Quals>
  </Read>
 </Chr>", "\n";
	}
	else { next; }
}
print STDOUT "</Contig>", "\n";
