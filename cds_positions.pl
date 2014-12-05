#!/uaopt/perl/5.14.2/bin/perl

use Bio::SeqIO;

my $stream = Bio::SeqIO->new(-file => "$ARGV[0]", -format => 'GenBank');
print STDERR "accessing file\n";
while (my $seq = $stream->next_seq() ) {
	print STDERR "opened seq\n";
	while (my $fth = $seq->_read_FTHelper_GenBank) {
		my $p_tag = $fth->key;
		print STDERR $p_tag, "\n";
		if ($p_tag eq "CDS") {
			# processing FTH into SeqFeqfeature object
			my $seqfeat= $fth->_generic_seqfeature($seq, "GenBank");
			my $start = $seqfeat->start;
			my $end = $seqfeat->end;
			my $strand = $seqfeat->strand;
			if ($strand == 1) {
				for (my $i = $start+2; $i <= $end; $i+=3) {
					print $i;
				}
			}
			if ($strand == -1) {
				for (my $i = $start; $i <= $end; $i+=3) {
					print $i;
				}
			}
		}
	}
}
