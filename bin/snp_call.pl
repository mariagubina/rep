#!/usr/bin/env perl

my ($patient, $bam_file, $genome) = @ARGV;
open(my $IN, "-|", "samtools view -q 20 -b $bam_file | bcftools mpileup -d 10000 -f $genome -");
$vcf_file = "${patient}.vcf";

open(my $OUT, ">", $vcf_file);

print $OUT join("\t", "CHR", "POS", "DP", "QS1", "QS2", "REF", "ALLELE1", "ALLELE2"), "\n";

while (my $line = <$IN>) {
    next if $line =~ /DEL/i;
    my @fields = split(/\t/, $line);

    if ($fields[4] =~ /,/) {
        $fields[7] =~ /DP=(\d+)/;
        my $dp = $1;

        if ($dp > 20) {
            $fields[7] =~ /QS=([e\-\d\.,]+)/;
            my @qs = split(/,/, $1);
            my @alt_alleles = split(/,/, $fields[4]);

            my @lit;
            my @QS;
            if ($qs[0] > 0.1) { push(@lit, $fields[3]); push(@QS, $qs[0]); }
            if ($qs[1] > 0.1) { push(@lit, $alt_alleles[0]); push(@QS, $qs[1]); }
            if ($qs[2] > 0.1) { push(@lit, $alt_alleles[1]); push(@QS, $qs[2]); }
            if ($qs[3] > 0.1) { push(@lit, $alt_alleles[2]); push(@QS, $qs[3]); }

            if (scalar(@lit) == 2) { 
                print $OUT join("\t", $fields[0], $fields[1], $dp, @QS, $fields[3], @lit), "\n";
            }
        }
    }
}
close($IN);
close($OUT);
