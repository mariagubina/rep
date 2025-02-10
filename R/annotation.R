library(readr)

library(Homo.sapiens)
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(GenomicRanges )

library(SNPlocs.Hsapiens.dbSNP150.GRCh38)
library(AllelicImbalance)
library(BSgenome.Hsapiens.UCSC.hg38)


stat_nf <- read_delim("/media/leon/DISK2/diabet/stat_nf.csv",
delim = "\t", escape_double = FALSE,
col_names = FALSE, trim_ws = TRUE, skip = 1)
View(stat_nf)
stat_nf<-na.omit(stat_nf)
colnames(stat_nf)=c("file","sample","chr","pos","ref","DP_in_sample","DP_in_file_ref","DP_in_file_alt","allel1","QS1_in_sample",     "QS1_in_ref","QS1_in_alt","allel2","QS2_in_sample","QS2_in_ref","QS2_in_alt")
stat_nf$QS1_in_ref<-as.numeric(as.character(stat_nf$QS1_in_ref))
stat_nf$cov_all1=stat_nf$DP_in_file_ref*stat_nf$QS1_in_ref
stat_nf$QS2_in_alt<-as.numeric(as.character(stat_nf$QS2_in_alt))
stat_nf$cov_all2=stat_nf$DP_in_file_alt*stat_nf$QS2_in_alt
stat_nf$DP_in_file_alt<-as.numeric(as.character(stat_nf$DP_in_file_alt))
stat_nf$cov_all2=stat_nf$DP_in_file_alt*stat_nf$QS2_in_alt


snp<-unique(stat_nf[,c(3,4,5)])
snp$chr<-as.character(as.numeric(snp$chr))
snp<-na.omit(snp)
SNPrenge<-GRanges(seqnames =as.character(snp$chr), IRanges(snp$pos,width = 1),ref=snp$ref)
updatedGRanges<-getSnpIdFromLocation(SNPrenge, SNPlocs.Hsapiens.dbSNP150.GRCh38)
updatedGRanges<-data.frame(updatedGRanges,rsid=names(updatedGRanges))
updatedGRanges<-na.omit(updatedGRanges)

########SNP_MAF<- getBM(attributes=c('chr_name',"start",'refsnp_id','minor_allele','minor_allele_freq'),filters = c('chr_name',"start"), values =as.list(snp[,c(1,2)]) , mart = snpMart)

#stat_shot<-stat_nf[,c(1,3,4,5,9,13,17,18)]
stat_shot<-merge(stat_nf,updatedGRanges[,c(1,2,7)],by.x=c(3,4),by.y=c(1,2))


library(biomaRt)
snpMart = useEnsembl(biomart = "snps",dataset = "hsapiens_snp")
SNP_MAF<- getBM(attributes=c('refsnp_id','minor_allele','minor_allele_freq'),filters = 'snp_filter', values = unique(stat_shot$rsid), mart = snpMart)
SNP_MAF=SNP_MAF[0,]
for (i in seq(1, length(unique(stat_shot$rsid)),by=1000)){SNP_MAF<-rbind(SNP_MAF,getBM(attributes=c('refsnp_id','minor_allele','minor_allele_freq'),filters = 'snp_filter', values = unique(stat_shot$rsid)[i:(i+1000)], mart = snpMart))}

stat_shot<-merge(stat_shot,SNP_MAF,by.x=19,by.y=1)
metadata <- read_delim("/media/leon/DISK2/diabet/metadata.txt",
delim = "\t", escape_double = FALSE,
col_names = FALSE, trim_ws = TRUE)
View(metadata)
stat_shot<-merge(stat_shot,metadata,by.x=4,by.y=3)



SNP_RNA<-stat_shot[stat_shot$X2=="RNA-seq",]
SNPrenge<-GRanges(seqnames = paste0("chr",as.character(SNP_RNA$chr)), IRanges(SNP_RNA$pos,width = 2))


TxDb(Homo.sapiens) <- TxDb.Hsapiens.UCSC.hg38.knownGene
exon<-exons(Homo.sapiens,columns=c("EXONID", "TXNAME", "GENEID","SYMBOL"))
x<-findOverlaps(exon,SNPrenge)

SNP_RNA<-SNP_RNA[subjectHits(x),]
SNP_RNA$gene<-sapply(exon[queryHits(x)]$GENEID, [, 1, simplify = TRUE)
SNP_RNA$Symbol<-sapply(exon[queryHits(x)]$SYMBOL, [, 1, simplify = TRUE)
SNP_RNA_<-unique(SNP_RNA[,c(1,2,5,6,10,14,18,19,20,26,27)])
library(plyr)
SNP_RNA_$cov_min_all<-NA

SNP_RNA_$cov_min_all[SNP_RNA_$allel1==SNP_RNA_$minor_allele]<-SNP_RNA_$cov_all1[SNP_RNA_$allel1==SNP_RNA_$minor_allele]
SNP_RNA_$cov_min_all[SNP_RNA_$allel2==SNP_RNA_$minor_allele]<-SNP_RNA_$cov_all2[SNP_RNA_$allel2==SNP_RNA_$minor_allele]
SNP_RNA_$cov_maj_all<-NA
SNP_RNA_$cov_maj_all[SNP_RNA_$allel1==SNP_RNA_$minor_allele]<-SNP_RNA_$cov_all2[SNP_RNA_$allel1==SNP_RNA_$minor_allele]
SNP_RNA_$cov_maj_all[SNP_RNA_$allel2==SNP_RNA_$minor_allele]<-SNP_RNA_$cov_all1[SNP_RNA_$allel2==SNP_RNA_$minor_allele]


ddply(SNP_RNA_,.(file,ref,allel1,allel2,minor_allele,Symbol,gene),summarize,cov_min=sum(cov_min_all),cov_maj=sum(cov_maj_all))
dfGene<-ddply(SNP_RNA_,.(file,minor_allele,Symbol,gene),summarize,cov_min=sum(cov_min_all),cov_maj=sum(cov_maj_all))
dfGene<-na.omit(dfGene)
dfGene$p.value<-apply(dfGene[,c(5,6)] ,1, function(x){binom.test(c(round(x[1]),round(x[2])),0.5)$p.value})
dfGene$p.adj<-p.adjust(dfGene$p.value,method = "BH")
listGene<-unique(dfGene[dfGene$p.adj<0.1,]$gene)

## ChIP gene annotation



SNP_ChIP<-stat_shot[stat_shot$X2!="RNA-seq",]
SNPrenge<-GRanges(seqnames = paste0("chr",as.character(SNP_ChIP$chr)), IRanges(SNP_ChIP$pos,width = 2))


##add enchansers

human_permissive_enhancers_phase_1_and_2 <- read.delim("~/Rproject/Alz/human_permissive_enhancers_phase_1_and_2.bed", header=FALSE)
SNPrenge<-as(human_permissive_enhancers_phase_1_and_2$V4,"GRanges")
x<-findOverlaps(promot,SNPrenge)
human_permissive_enhancers_phase_1_and_2<-human_permissive_enhancers_phase_1_and_2[subjectHits(x),]
human_permissive_enhancers_phase_1_and_2$GENEID<-as.character(promot[queryHits(x)]$GENEID)
human_permissive_enhancers_phase_1_and_2$Symbol<-as.character(promot[queryHits(x)]$SYMBOL)
human_permissive_enhancers_phase_1_and_2$TXNAME<-as.character(promot[queryHits(x)]$TXNAME)
human_permissive_enhancers_phase_1_and_2<-unique(human_permissive_enhancers_phase_1_and_2)
human_permissive_enhancers_phase_1_and_2<-na.omit(human_permissive_enhancers_phase_1_and_2)


enchans<-GRanges(seqnames =as.character(human_permissive_enhancers_phase_1_and_2$V1), IRanges(start = human_permissive_enhancers_phase_1_and_2$V2,end = human_permissive_enhancers_phase_1_and_2$V3),GENEID=human_permissive_enhancers_phase_1_and_2$GENEID,SYMBOL=human_permissive_enhancers_phase_1_and_2$Symbol,TXNAME=human_permissive_enhancers_phase_1_and_2$TXNAME,type="enchancer")






TxDb(Homo.sapiens) <- TxDb.Hsapiens.UCSC.hg38.knownGene
gene<-transcripts(Homo.sapiens,columns=c("EXONID", "TXNAME", "GENEID","SYMBOL"))
promot<-c(GRanges(seqnames = seqnames(gene), IRanges(start=start(gene)-1000,width = 2000),TXNAME=gene$TXNAME,SYMBOL=gene$SYMBOL,GENEID=gene$GENEID),GRanges(seqnames = seqnames(gene), IRanges(start=end(gene)-1000,width = 2000),GENEID=gene$GENEID,SYMBOL=gene$SYMBOL,TXNAME=gene$TXNAME))
promot$type="promoter"
promot<-c(promot,enchans)


x<-findOverlaps(promot,SNPrenge)

SNP_ChIP<-SNP_ChIP[subjectHits(x),]
SNP_ChIP$gene<-as.character(promot[queryHits(x)]$GENEID)
SNP_ChIP$Symbol<-as.character(promot[queryHits(x)]$SYMBOL)
SNP_ChIP$regtype<-as.character(promot[queryHits(x)]$type)
SNP_ChIP<-unique(SNP_ChIP[,c(1,2,5,6,10,14,18,19,20,26,27)])
SNP_ChIP$chr<-as.numeric(SNP_ChIP$chr)
SNP_ChIP<-na.omit(SNP_ChIP)
SNP_ChIP$chr<-as.character(SNP_ChIP$chr)
SNPrenge<-GRanges(seqnames =as.character(SNP_ChIP$chr), IRanges(SNP_ChIP$pos,width = 1),ref=SNP_ChIP$reference)
updatedGRanges<-getSnpIdFromLocation(SNPrenge, SNPlocs.Hsapiens.dbSNP151.GRCh38)
updatedGRanges<- data.frame(chr=seqnames(updatedGRanges),pos=start(updatedGRanges),snpid=names(updatedGRanges))
SNP_ChIP<-merge(SNP_ChIP,unique(updatedGRanges),by.x=c(2,3),by.y=c(1,2))
SNP_ChIP<-na.omit(SNP_ChIP)
SNP_ChIP<-SNP_ChIP[is.element(SNP_ChIP$gene,listGene),]

SNP_ChIP$p.value<-apply(SNP_ChIP[,c(7,8)] ,1, function(x){binom.test(c(round(x[1]),round(x[2])),0.5)$p.value})
SNP_ChIP$p.adj<-p.adjust(SNP_ChIP$p.value,method = "BH")
SNP_reg<-unique(SNP_ChIP[SNP_ChIP$p.adj<0.1,])

GWAS <- read.delim("~/Rproject/Alz/GWAS")
SNP_ChIP_<-SNP_ChIP[SNP_ChIP$p.adj<0.1,]
SNP_reg<-unique(SNP_ChIP[SNP_ChIP$p.adj<0.1,])


TXN <- read_delim("~/Rproject/Alz/TXN", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
TXN<-na.omit(TXN)
SNP_reg_pos<-unique(stat_shot[is.element(stat_shot$rsid,unique(SNP_reg$rsid)),c(2,3,4)])


TXNrenge<-GRanges(seqnames =as.character(TXN$chrom), IRanges(start=as.numeric(TXN$chromStart),end = TXN$chromEnd),Symbol=TXN$name)
SNPrenge<-GRanges(seqnames =paste0("chr",as.character(SNP_reg_pos$chr)), IRanges(start=as.numeric(SNP_reg_pos$pos),end = SNP_reg_pos$pos+2),Symbol=SNP_reg_pos$rsid)

x<-findOverlaps(SNPrenge,TXNrenge)
SNP_reg_pos<-SNP_reg_pos[queryHits(x),]
SNP_reg_pos$TF<-TXNrenge[subjectHits(x)]$Symbol
SNP_reg_pos<-unique(SNP_reg_pos)
SNP_TXN_reg<-(unique(SNP_reg_pos$rsid))

library(motifbreakR)

library(SNPlocs.Hsapiens.dbSNP150.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38)
x<-SNP_TXN_reg
k=0
for (i in unique(x)){
  k=k+1
  variant=snps.from.rsid(rsid = i,
                         dbSNP = SNPlocs.Hsapiens.dbSNP150.GRCh38,
                         search.genome = BSgenome.Hsapiens.UCSC.hg38)
  result <- motifbreakR(snpList = variant,pwmList = MotifDb, threshold = 0.9)
 all.result[k]<-result
try(names(all.result[k])<-i)
 }
SNP_reg<-merge(SNP_reg,SNP_reg_pos,by="rsid",all.x=T)
