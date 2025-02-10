#!/usr/bin/env perl

my (@stats) = @ARGV;

open(my $OUT_ref, ">", "ref.stat");
open(my $OUT_alt, ">", "alt.stat");

print $OUT_ref join("\t", "FILE", "SAMPLE", "CHR", "POS", "DP_ref", "QS1_ref", "QS2_ref", "REF", "ALLELE1_ref", "ALLELE2_ref"), "\n";
print $OUT_alt join("\t", "FILE", "SAMPLE", "CHR", "POS", "DP_alt", "QS1_alt", "QS2_alt", "REF", "ALLELE1_alt", "ALLELE2_alt"), "\n";

foreach my $file (@stats) {
    open(my $IN, "<", $file);
    <$IN>;  # Skip header

    if ($file =~ /alt/) {
        while (my $line = <$IN>) { print $OUT_alt $line; }
    } else {
        while (my $line = <$IN>) { print $OUT_ref $line; }
    }
    close($IN);
}
close($OUT_ref);
close($OUT_alt);
