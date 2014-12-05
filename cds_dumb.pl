my $site=$ARGV[0];

while (<STDIN>) {
	if ($_ =~ 'CDS') {
		if ($_ =~ 'complement') {
			$_ =~ m/.*?(\d+)\.\.(\d+)/;
			#print "complement ", $1, "..", $2, "\n";
			for (my $i = $1 + 3 - $site; $i <= $2; $i+=3) {
				print $i, "\n";
			}
		}
		else {
			$_ =~ m/.*?(\d+)\.\.(\d+)/;
			#print $1, "..", $2, "\n";
			for (my $i = $1 + $site - 1; $i <= $2; $i+=3) {
				print $i, "\n";
			}
		}
	}
}
