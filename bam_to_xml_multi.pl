#!/uaopt/perl/5.14.2/bin/perl

# A quick converter of aligned reads in sam format to xml, for use with PIIM
# Assumes the input bam file has been sorted; otherwise may take extra time to sort while processing
# Based largely on examples provided by Lincoln Stein at http://search.cpan.org/~lds/Bio-SamTools/lib/Bio/DB/Sam.pm

use warnings;
use strict;
use Bio::DB::Sam;
use Bio::SeqIO;
use List::Util qw(sum);

# a very short subroutine to calculate the mean of elements in an array
sub mean {
	return sum(@_)/@_;
}
  
# load bam file specified at the command line

my $sam = Bio::DB::Sam->new(-bam  =>"$ARGV[0]",
                             -fasta=>"$ARGV[1]");
my @targets = $sam->seq_ids;
							 
GENOME:foreach my $target (@targets) {
	if (-e "$target.xml") {
		print "File $target.xml exists; skipping $target reads.\n";
		next GENOME;
	}
	open (OUTFILE, ">$target.xml");
	# print header information to outfile
	print OUTFILE "<?xml version=\"1.0\"?>
<!DOCTYPE Contig
SYSTEM \"piim_r.dtd\">
<Contig>", "\n";

	# error file specifies options used to run
	open (ERRFILE, ">$target.err");
	
	# a little note to the error file on the provenance of the data
	print ERRFILE "Processing the pairs of reads from $ARGV[0] which mapped to genomes in $ARGV[1].\n";
	if ($ARGV[3] eq 'min') {
		print ERRFILE "Filtering down to minimum $ARGV[2] bp coverage.\n";
	} elsif ($ARGV[3] eq 'max') {
		print ERRFILE "Filtering down to maximum $ARGV[2] bp coverage.\n";
	} else { 
		die "Usage: bam_to_matrix.pl [in.bam] [in.fasta] [coverage] [min/max]\n";
	}
	
	# load all pairs into a single array; may be the memory-limiting step
	my @pairs = $sam->get_features_by_location(-type   => 'read_pair',
						-seq_id => "$target");
	
	# introduce an array to keep track of coverage levels across the genome
	my @coverage;
		
	# update .xml file one pair of reads at a time
	PAIRLOOP: for my $pair (@pairs) {
		
		# grab the pair of reads, check that both map and get locations on the reference
		my ($first_mate,$second_mate) = $pair->get_SeqFeatures;
		if ($first_mate->proper_pair) {
			my $f_start = $first_mate->start;
			my $s_start = $second_mate->start;
			my $f_end = $first_mate->end;
			my $s_end = $second_mate->end;
			
			# check if outputting the pair of reads would exceed specified coverage limit
			if ($ARGV[3] eq 'max') {
				foreach my $i ($f_start..$f_end) {
					if ($coverage[$i]>$ARGV[2]) {
						next PAIRLOOP; # proceed to next pair
					}
				}
				foreach my $i ($s_start..$s_end) {
					if ($coverage[$i]>$ARGV[2]) {
						next PAIRLOOP;
					}
				}
			}
			elsif ($ARGV[3] eq 'min') {
				my $f_ok = 0;
				foreach my $i ($f_start..$f_end) {
					if ($coverage[$i]<$ARGV[2]) {
						$f_ok = 1;
					}
				}
				unless ($f_ok) {
					next PAIRLOOP;
				}
				my $s_ok = 0;
				foreach my $i ($s_start..$s_end) {
					if ($coverage[$i]<$ARGV[2]) {
						$s_ok = 1;
					}
				}
				unless ($s_ok) {
					next PAIRLOOP;
				}
			}

			# update genome coverage once both reads are accepted for output
			foreach my $i ($f_start..$f_end) {
				$coverage[$i]++;
			}
			foreach my $i ($s_start..$s_end) {
				$coverage[$i]++;
			}
			
			# gather read data
			my $f_name = $first_mate->name;
			my $s_name = $second_mate->name;
			my $f_dna = $first_mate->query->dna;
			my $s_dna = $second_mate->query->dna;
			my @f_scores = $first_mate->qscore;
			my @s_scores = $second_mate->qscore;
			
			# print header info for first read
			print OUTFILE " <Chr>
  <Read start=\"$f_start\">
   <Bases>";

			# insert length may be an indicator of mapping quality
			my $length = $pair->length;
			print ERRFILE "read pair $f_name and $s_name located between $f_start and $s_end with insert length $length\n";
			
			# reading CIGAR string to exclude indels from nucleotide count matrix
			my $cigar = $first_mate->cigar_str;
			my $cigar_len=0;
			my $cigar_offset=0;
			my @match_quals;
			
			# process CIGAR string one piece at a time, from start to finish
			while ($cigar=~s/(\d+\D)(.*)/$2/) {
				my $cig=$1;
				$cig=~s/(\d+)(\D)/$2/;
				my $cignum=$1;
				print ERRFILE "processing first read $cig region $cignum bp long, starting at ", $f_start+$cigar_offset, "\n";
				if ($cig eq 'M') {
					foreach my $i (0..$cignum-1) {
						if (substr($f_dna,$cigar_len+$i,1) eq 'A') {
							print OUTFILE "A";
							print ERRFILE "A";
						}
						if (substr($f_dna,$cigar_len+$i,1) eq 'T') {
							print OUTFILE "T";
							print ERRFILE "T";
						}
						if (substr($f_dna,$cigar_len+$i,1) eq 'C') {
							print OUTFILE "C";
							print ERRFILE "C";
						}
						if (substr($f_dna,$cigar_len+$i,1) eq 'G') {
							print OUTFILE "G";
							print ERRFILE "G";
						}
					}
					foreach my $i ($cigar_len..$cigar_len+$cignum-1) {
						push @match_quals, $f_scores[$i]; 
					}
					$cigar_len+=$cignum;
					$cigar_offset+=$cignum;
					print ERRFILE "\nprocessed $cignum bp and now at read position $cigar_len\n";
				}
				elsif ($cig eq 'I') {
					my @insert_quals;
					foreach my $i ($cigar_len..$cigar_len+$cignum-1) {
						push @insert_quals, $f_scores[$i];
						print ERRFILE substr($f_dna, $i, 1); 
					}
					my $insert_qual = mean(@insert_quals);
					print ERRFILE "\ninsertion in first read $f_name is $cigar_len bp from the start, at ", $f_start+$cigar_offset, " bp, ",scalar(@insert_quals)," bp long and average quality of $insert_qual\n";
					$cigar_len+=$cignum;
				}
				elsif ($cig eq 'S') {
					foreach my $i ($cigar_len..$cigar_len+$cignum-1) {
						print ERRFILE substr($s_dna, $i, 1);
					}
					print ERRFILE "\nsoft clipped in second read $s_name is $cigar_len bp from the start, at ", $s_start+$cigar_offset, " bp, ",$cignum," bp long\n";
					$cigar_len+=$cignum;
				}
				elsif ($cig eq 'D') {
					print ERRFILE "deletion in first read $f_name is $cigar_len bp from the start, at ", $f_start+$cigar_offset, " bp, $cignum bp long\n";
					$f_start+=$cignum;
				}
				else { print ERRFILE "unexpected CIGAR string for first read $f_name: $cignum, $cig\n"; }
			}
			if ($cigar ne '' || $cigar_len != scalar(@f_scores)) { 
				print ERRFILE "CIGAR string mishandled for first read $f_name: $cigar_len bp processed out of ", scalar(@f_scores), "\n";
			}
			
			# finishing output of first read and starting the second
			print OUTFILE "</Bases>
   <Quals>@match_quals</Quals>
  </Read>
  <Read start=\"$s_start\">
   <Bases>";
   
			# reading CIGAR string to exclude indels from nucleotide count matrix
			$cigar = $second_mate->cigar_str;
			$cigar_len=0;
			$cigar_offset=0;
			my @match_quals;
			
			# process CIGAR string one piece at a time, from start to finish
			while ($cigar=~s/(\d+\D)(.*)/$2/) {
				my $cig=$1;
				$cig=~s/(\d+)(\D)/$2/;
				my $cignum=$1;
				print ERRFILE "processing second read $cig region $cignum bp long, starting at ", $s_start+$cigar_offset, "\n";
				if ($cig eq 'M') {
					foreach my $i (0..$cignum-1) {
						if (substr($s_dna,$cigar_len+$i,1) eq 'A') {
							print OUTFILE "A";
							print ERRFILE "A";
						}
						if (substr($s_dna,$cigar_len+$i,1) eq 'T') {
							print OUTFILE "T";
							print ERRFILE "T";
						}
						if (substr($s_dna,$cigar_len+$i,1) eq 'C') {
							print OUTFILE "C";
							print ERRFILE "C";
						}
						if (substr($s_dna,$cigar_len+$i,1) eq 'G') {
							print OUTFILE "G";
							print ERRFILE "G";
						}
					}
					foreach my $i ($cigar_len..$cigar_len+$cignum-1) {
						push @match_quals, $f_scores[$i]; 
					}
					$cigar_len+=$cignum;
					$cigar_offset+=$cignum;
					print ERRFILE "\nprocessed $cignum bp and now at read position $cigar_len\n";
				}
				elsif ($cig eq 'I') {
					my @insert_quals;
					foreach my $i ($cigar_len..$cigar_len+$cignum-1) {
						push @insert_quals, $s_scores[$i];
						print ERRFILE substr($s_dna, $i, 1);
					}
					my $insert_qual = mean(@insert_quals);
					print ERRFILE "\ninsertion in second read $s_name is $cigar_len bp from the start, at ", $s_start+$cigar_offset, " bp, ",scalar(@insert_quals)," bp long and average quality of $insert_qual\n";
					$cigar_len+=$cignum;
				}
				elsif ($cig eq 'D') {
					print ERRFILE "deletion in second read $s_name is $cigar_len bp from the start, at ", $s_start+$cigar_offset, " bp, $cignum bp long\n";
					$s_start+=$cignum;
				}
				elsif ($cig eq 'S') {
					foreach my $i ($cigar_len..$cigar_len+$cignum-1) {
						print ERRFILE substr($s_dna, $i, 1);
					}
					print ERRFILE "\nsoft clipped in second read $s_name is $cigar_len bp from the start, at ", $s_start+$cigar_offset, " bp, ",$cignum," bp long\n";
					$cigar_len+=$cignum;
				}
				else { print ERRFILE "unexpected CIGAR string for second read $s_name: $cignum, $cig\n"; }
			}
			if ($cigar ne '' || $cigar_len != scalar(@s_scores)) { 
				print ERRFILE "CIGAR string mishandled for second read $s_name: $cigar_len bp processed out of ", scalar(@s_scores), "\n";
			}
			
			# finish output of second read
			print OUTFILE "</Bases>
   <Quals>@match_quals</Quals>
  </Read>
 </Chr>", "\n";
		}
		else { next; }
	}
	print OUTFILE "</Contig>", "\n";
}
