#!/usr/bin/perl
($sample, $bam_file, $genome) = @ARGV;
$sample =~ s/^[^_]*_//; 
open(OUT,">$sample.stat");
print OUT join("\t", "SAMPLE", "CHR", "POS", "DP_total", "REF", "QS_ref", "ALT1", "QS_alt1", "ALT2", "QS_alt2"), "\n";
# open(IN,"samtools view -bq 30 $samp |bcftools mpileup -f /media/leon/DATA/Genomes/grch38_snp/Human_hg38.fa -d 10000 -|");
# $genome = "/media/leon/Polina/Genomes/hg38.chromFa/hg38_FOR_HISAT_OUTPUT.fa"

open(IN,"samtools view -bq 30 $bam_file |bcftools mpileup -f $genome -d 10000 -|");

while ($IN=<IN>) {
#	print $IN;
	$IN=~/DP=(\d+)/;
	$dp=$1;
	if ($dp>50) {
		@IN=split(/\t/,$IN);
		if($IN[4]=~/,/) {
			$IN[7]=~/QS=([\d\.,]+)/;
			@QS=split(/,/,$1);
			if($QS[1]>0.1 && $QS[1]<0.9) {
				@l=split(/,/,$IN[4]);
				print OUT	"$sample	$IN[0]	$IN[1]	$dp	$IN[3]	$QS[0]	$l[0]	$QS[1]	$l[1]	$QS[2]\n";
			}
		} 
	}
}
