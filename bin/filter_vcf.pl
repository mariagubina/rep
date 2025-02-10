#!/usr/bin/env perl

my ($patient, $sample, $bam_file, $vcf_file) = @ARGV;

my $ver = "ref";
if ($bam_file =~ /alt/) { $ver = "alt"; }

open(my $OUT, ">", "${patient}_${sample}_${ver}.stat"); # создаём файл для записи результатов для данного bam-файла (ref или alt для определённого образца)
# что значит _ref? (должно быть _$ver ?)
print $OUT join("\t", "FILE", "SAMPLE", "CHR", "POS", "DP_ref", "QS1_ref", "QS2_ref", "REF", "ALLELE1_ref", "ALLELE2_ref"), "\n";

my %ref; # записываем весь vcf-файл для этого человека в хэш: хромосома+позиция =>  целая строка 
open(my $IN, "<", $vcf_file) or die "Failed to open VCF file $vcf_file: $!";
while (my $IN = <$IN>) {
    chomp $IN;
    my @IN = split(/\t/, $IN);
    $ref{"$IN[0]\t$IN[1]"} = $IN;
}
close($IN);

# вызываем пайлап для данного бам-файла
open(my $vcf_stream, "-|", "samtools view -q20 -b $bam_file | bcftools mpileup -d10000 --no-reference -");
while (my $line = <$vcf_stream>) {
    chomp $line;
    my @IN = split(/\t/, $line); # tab-delimitted строка из vcf в виде массива
    # для каждой позиции из пайлапа если она представлена в общем всф файле для данного человека
    # разбираем поля и записываем в $sample.stat_$gm.csv
    if ($ref{"$IN[0]\t$IN[1]"} ne "") { 
        $IN[7] =~ /DP=(\d+)/; # извлекаем глубину покрытия
        my $dp = $1;

        $IN[7] =~ /QS=([\d.,]+)/; # извлекаем строку с кволити скорами
        my @qs = split(/,/, $1); # делим по запятым и записываем в массив

        $IN[4] =~ s/\<\*\>/$IN[3]/; # если нет альтернативного аллеля то берём референсный?
        my @alt = split(/,/, $IN[4]); # если есть альт аллели делим по запятым и записываем в массив

        my @lit; # массив со всеми аллелями
        my @QS; # массив с кволити скорами для всех аллелей
        push(@lit, $IN[3]), push(@QS, $qs[0]);
        push(@lit, $alt[0]), push(@QS, $qs[1]);
        push(@lit, $alt[1]), push(@QS, $qs[2]);
        push(@lit, $alt[2]), push(@QS, $qs[3]);

        print $OUT join("\t", $sample, $patient, $IN[0], $IN[1], $dp, @QS, $IN[3], @lit), "\n";
    }
}
close($vcf_stream);
close($OUT);
