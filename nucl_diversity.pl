#!/uaopt/perl/5.14.2/bin/perl

# Quick script for calculating pi and counting the number of variable sites
#
# Input: nucleotide count matrix (from specified file)
# Output: number of variable sites (to std error), pi per site (to std output)

use List::Util qw(max);

my $var_sites;
my $tot_sites;
my $genome_len;
my $div_tot;

open (INF, $ARGV[0]);

while (<INF>) {
	$genome_len++;
	chomp $_;
	my ($pos,$a,$t,$c,$g,$cov)=split('\t', $_);
	if ($cov > 0) {
		$tot_sites++;
		my $a_frq = $a/$cov;
		my $t_frq = $t/$cov;
		my $c_frq = $c/$cov;
		my $g_frq = $g/$cov;
		my $maj_frq = max $a_frq,$t_frq,$c_frq,$g_frq;
		# making a folded allele frequency spectrum
		$maj_frq = 1 - $maj_frq;
		# in case major allele frequency is below 0.5
		if ($maj_frq > 0.5) { $maj_frq = 1 - $maj_frq; }
		my $diversity = $a_frq*($c_frq+$g_frq+$t_frq) + $c_frq*($g_frq+$t_frq) + $g_frq*$t_frq;
		$div_tot += $diversity;
		print STDOUT "$pos", "\t", "$diversity", "\t", "$cov", "\t",  "$maj_frq", "\n";
		if ($diversity != 0) {
			$var_sites++;
		}
	}
}
my $div_avg = $div_tot / $tot_sites;
print STDERR "Analyzed $tot_sites sites and found $var_sites variable sites and average nucleotide diversity of $div_avg for a $genome_len bp genome of ", $ARGV[0], "\n";
