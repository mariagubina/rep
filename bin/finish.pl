#!/usr/bin/env perl

my ($stat_ref, $stat_alt, @vcfs) = @ARGV;

open(my $OUT, ">", "stat_nf.csv");
print $OUT join("\t", "FILE", "SAMPLE", "CHR", "POS", "REF", "DP_in_sample", "DP_in_file",
                "ALLELE1", "QS1_in_sample", "QS1_in_file_ref", "QS1_in_file_alt",
                "ALLELE2", "QS2_in_sample", "QS2_in_file_ref", "QS2_in_file_alt"), "\n";

my %ref; # хэш пациент+хромосома+позиция => вся строка (один для всех 34 всф файлов)

foreach my $file (@vcfs) { # перебираем каждый VCF файл
    my ($patient) = $file =~ m|([^/]+)\.vcf$|;
    open(my $VCF_IN, "<", $file);
    # print "$patient\n$file\n"; # можно добавить тэг в директиве процесса в некстфлоу
    while (my $line = <$VCF_IN>) {
        chomp $line;
        my @fields = split(/\t/, $line);
        my $key = "$patient\t$fields[0]\t$fields[1]";
        $ref{$key} = $line;
    }
}

open(my $ref_in, "<", $stat_ref); # открываем смёрдженную таблицу со статистикой по снп для всех бам-файлов (реф)
while (my $IN = <$ref_in>) {
    chomp($IN);
    @IN = split(/\t/, $IN);
    my $key = "$IN[1]\t$IN[2]\t$IN[3]"; # "SAMPLE", "CHR", "POS"
    if (exists $ref{$key}) {
        chomp($ref{$key});
        @vcf = split(/\t/, $ref{$key}); # "CHR", "POS", "DP", "QS1", "QS2", "REF", "ALLELE1", "ALLELE2"
        $a1{$IN[11]} = $IN[6]; # ALT1 => QS1
        $a1{$IN[12]} = $IN[7]; # ALT2 => QS2
    
        my $key = join("\t", "$IN[0]", "$IN[1]", "$IN[2]", "$IN[3]");
        #                      FILE     SAMPLE     CHR       POS
        my $value = join("\t", "$vcf[5]", "$vcf[2]", "$IN[4]", "$vcf[6]", "$vcf[3]", "$a1{$vcf[6]}", "$vcf[7]", "$vcf[4]", "$a1{$vcf[7]}");
        #                        REF   DP_in_sample DP_in_file ALLELE1  QS1_in_sample   QS1_in_file  ALLELE2 QS2_in_sample QS2_in_file   
        $stat{$key} = $value; # записывем в хэш stat по ключу файл+сэмпл+хромосома+позиция метрики из референсного бам файла
        undef %a1;
    }
}

open(my $alt_in, "<", $stat_alt);

while (my $IN = <$alt_in>) {
    chomp($IN);
    @IN = split(/\t/, $IN);
    my $key = "$IN[1]\t$IN[2]\t$IN[3]";
    if (exists $ref{$key}) {
        $a1{$IN[11]} = $IN[6];
        $a1{$IN[12]} = $IN[7];

        my $key = join("\t", "$IN[0]", "$IN[1]", "$IN[2]", "$IN[3]");
        @refs = split(/\t/, $stat{$key}); # находим в хэше stat метрики для данной позиции для референсного бам файла
        # переписываем значение: REF, DP_in_sample, DP_in_file, DP_in_file опять, ALLELE1, QS1_in_sample, 
        my $value = join("\t", "$refs[0]", "$refs[1]", "$refs[2]", "$IN[4]", "$refs[3]",
        
                            "$refs[4]", "$refs[5]", "$a1{$refs[3]}", "$refs[6]", "$refs[7]",
                            "$refs[8]", "$a1{$refs[6]}");
        $stat{$key} = $value;
        undef %a1;
    }
}
foreach $key (keys(%stat)) { print $OUT "$key\t$stat{$key}\n" };
