#!/usr/bin/env perl

my ($patient, $vcf_file, $ref_genome) = @ARGV;
my (%GN, @chr, $chr, @seq);

open(my $IN, "<", $ref_genome);
while (my $line = <$IN>) {
    chomp($line);
    if ($line =~ /^>/) {
        if ($chr) {
            $GN{$chr} = join("", @seq);
            undef @seq;
        }
        $line =~ /(>chr\S+)/;
        $chr = $1;
        push(@chr, $chr);
    } else {
      push(@seq, $line);
    }
}
$GN{$chr} = join("", @seq) if $chr;
close($IN);

open(my $IN, "<", $vcf_file);
<$IN>;
while (my $line = <$IN>) {
    chomp($line);
    my @VC = split(/\t/, $line);
    if ($VC[5] eq $VC[6]) {
        my $rf = substr($GN{">$VC[0]"}, $VC[1] - 1, 1);
        substr($GN{">$VC[0]"}, $VC[1] - 1, 1) = $VC[7];
            
        if ($rf ne $VC[5]) {
            print "$rf\t$VC[5]\t$VC[7]\n";
        }
    }
}
close($IN);

$genome_file = "${patient}/${patient}.fa";
open(my $out, ">", $genome_file);
foreach my $key (@chr) {
    print $out "$key\n$GN{$key}\n";
}
close($out);
