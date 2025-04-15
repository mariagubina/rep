#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.bowtie_fa     = '/media/leon/Polina/Genomes/hg38.chromFa/hg38.fa'
params.bowtie_index  = '/media/leon/Polina/Genomes/hg38.chromFa/hg38'
params.hisat_fa      = '/media/leon/Polina/Genomes/hg38.chromFa/hg38_FOR_HISAT.fa'
params.hisat_index   = '/media/leon/DISK2/icig/grch38_snp_tran/genome_snp_tran'
params.metadata      = 'metadata.csv'
params.reads         = 'data/*{1,2}.fq.gz'
params.atac_dir      = '/media/leon/Masha/ATAC/fastqs'
params.mm_dir        = '/media/leon/Polina/atac_rna/raw_fastqs'
params.barcodes_dir  = '/media/leon/Masha/barcodes'
params.outdir        = 'results'
params.threads       = 4
params.more_threads  = 6
params.atac_adapters = '/usr/share/trimmomatic/NexteraPE-PE.fa'
params.gtf           = '/media/leon/DISK2/icig/Homo_sapiens.GRCh38.113.gtf'


process FASTQC {
    publishDir "${params.outdir}/QC", mode: 'copy'

    input:
    tuple val(sample), path(fastqs)

    output:
    tuple val(sample), path("*_fastqc.html"), path("*_fastqc.zip")

    script:
    """
    fastqc ${fastqs}
    """
}

process PREPROCESS {
    maxForks 8
    cpus 4
    memory '4GB'
    publishDir "${params.outdir}/trimmed", pattern: "*.fastq.gz"
    publishDir "${params.outdir}/trimmed/fastp_reports", pattern: "*.{html,json}"
    
    input:
    tuple val(sample), path(fastqs), val(assay_type)
    
    output:
    tuple val(sample), path("*.trimmed.fastq.gz"), val(assay_type), emit: reads
    tuple path("*.json"), path("*.html"), emit: reports
    
    script:
    if (assay_type == "ATAC" || assay_type == "mmATAC") {
        """
        umi_skip=""
        if [ ${assay_type} == "mmATAC" ]; then
            umi_skip="--umi_skip 8"
        fi

        ## Attaching barcodes to forward read headers
        fastp -i ${fastqs[0]} -I ${fastqs[1]} \\
            -o ${sample}_1.barcoded.fastq.gz -O tmp.fq.gz \\
            --umi --umi_loc=read2 ${umi_skip} --umi_len=16 --umi_prefix=CR \\
            --disable_length_filtering --disable_adapter_trimming --disable_quality_filtering \\
            --thread ${task.cpus}
        
        ## Attaching barcodes to reverse read headers
        fastp -i ${fastqs[2]} -I ${fastqs[1]} \\
            -o ${sample}_2.barcoded.fastq.gz -O tmp.fq.gz \\
            --umi --umi_loc=read2 ${umi_skip} --umi_len=16 --umi_prefix=CR \\
            --disable_length_filtering --disable_adapter_trimming --disable_quality_filtering \\
            --thread ${task.cpus}
        
        ## Trimming adapters and filtering
        fastp -i ${sample}_1.barcoded.fastq.gz -I ${sample}_2.barcoded.fastq.gz \\
            -o ${sample}_1.trimmed.fastq.gz -O ${sample}_2.trimmed.fastq.gz \\
            --length_required 36 --trim_poly_g \\
            --adapter_fasta ${params.atac_adapters} \\
            --thread ${task.cpus} \\
            -j ${sample}.fastp.json -h ${sample}.fastp.html
        
        ## Removing intermediate files
        rm tmp.fq.gz ${sample}_1.barcoded.fastq.gz ${sample}_2.barcoded.fastq.gz
        """
    } else if (assay_type == "mmRNA") {
        """
        fastp \\
            -i ${fastqs[0]} -I ${fastqs[1]} \\
            -o ${sample}_1.trimmed.fastq.gz -O ${sample}.trimmed.fastq.gz \\
            --umi --umi_loc=read1 --umi_len=16 --umi_prefix=CR \\
            --length_required 36 --trim_poly_g \\
            -a AAGCAGTGGTATCAACGCAGAGTAC \\ # сделать mmRNA_adapter параметром
            --thread ${task.cpus} \\
            -j ${sample}.json -h ${sample}.html  
        
        rm ${sample}_1.trimmed.fastq.gz
        """
    } else {
        """
        fastp \\
            -i ${fastqs[0]} -I ${fastqs[1]} \\
            -o ${sample}_1.trimmed.fastq.gz -O ${sample}_2.trimmed.fastq.gz \\
            --detect_adapter_for_pe --trim_poly_g \\
            --thread ${task.cpus} \\
            --length_required 36 \\
            -j ${sample}.fastp.json -h ${sample}.fastp.html
        """
    }
}

process HISAT2 {
    cpus 30
    memory '30GB'
    maxForks 1    
    publishDir "${params.outdir}/alignments"

    input:
    tuple val(sample), path(reads), val(assay_type)

    output:
    tuple val(sample), path("${sample}.bam"), val(assay_type), emit: bam
    tuple val(sample), path("${sample}.bam.bai"), val(assay_type), emit: bai

    script:
    """
    splicing_mode="--no-spliced-alignment" # пусть лучше условие будет вне кавычек (?)
    if [ ${assay_type} == "RNASEQ" ]; then
        splicing_mode=""
    fi

    hisat2 -p ${task.cpus} -x ${params.hisat_index} \\
        -1 ${reads[0]} -2 ${reads[1]} \${splicing_mode} \\
        | samtools view -@ ${task.cpus} -b -u - \\
        | samtools sort -@ ${task.cpus} -o ${sample}.bam

    samtools index -@ ${task.cpus} ${sample}.bam
    """
}

process FEATURE_COUNTS {
    cpus 20
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    path(alignments)

    output:
    tuple path("counts.tsv"), path("counts.tsv.summary")

    script:
    """
    featureCounts -a ${params.gtf} -o counts.tsv -T ${task.cpus} -p -B ${alignments}
    """
}

process GET_STATS {
    cpus 1
    memory '3GB'
    maxForks 10 

    publishDir "${params.outdir}/stats", mode: 'copy'

    input:
    tuple val(sample), path(bam_file), val(assay_type)

    output:
    tuple val(sample), path("*.stat")

    script:
    """
    stats.pl ${sample} ${bam_file} ${params.hisat_fa}
    """
}

// process ANNOTATE {
//     publishDir "${params.outdir}/stats", mode: 'move'

//     input:
//     tuple path(stats_list)

//     output:
//     path("*.csv")

//     script:
//     """
//     asym.R ${params.metadata} ${stats_list}
//     """
// }

process BOWTIE2 {
    cpus 15
    memory '20GB'
    maxForks 2

    publishDir "${params.outdir}/ref"

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}.bam"), emit: bam
    tuple val(sample), path("${sample}.bam.bai"), emit: bai

    script:
    def bowtie_cmd = reads.size() == 1 ? 
        "-U ${reads[0]}" :
        "-1 ${reads[0]} -2 ${reads[1]}"
    """
    bowtie2 --no-unal -p ${task.cpus} -x ${params.bowtie_index} ${bowtie_cmd} \\
        | samtools view -@ ${task.cpus} -b -q 10 -u - \\
        | samtools sort -@ ${task.cpus} -o ${sample}.bam

    samtools index -@ ${task.cpus} ${sample}.bam
    """
}

process MERGE_BAMS {
    publishDir "${params.outdir}/merged"

    input:
    tuple val(patient), path(bams)

    output:
    tuple val(patient), path("${patient}.bam"), emit: bam
    path("${patient}.bam.bai"), emit: bai

    script:
    """
    if [[ " MM120 MM121 MM124 MM54 MM56 MM57 MM79 MM80 MM81 MM86 MM94 JYH792 JYH809 MM12 " =~ " ${patient} " ]]; then
        ln -s ${bams} ${patient}.bam
        samtools index -@ ${task.cpus} ${patient}.bam
    else
        samtools merge -@ ${task.cpus} ${patient}.bam ${bams}
        samtools index -@ ${task.cpus} ${patient}.bam
    fi
    """
}

process SNP_CALL {
    publishDir "${params.outdir}/merged"

    input:
    tuple val(patient), path(merged_bam)

    output:
    tuple val(patient), path("${patient}.vcf")

    script:
    """
    snp_call.pl ${patient} ${merged_bam} ${params.bowtie_fa}
    """
}

process ALT_GENOME {
    publishDir "${params.outdir}/genomes"

    input:
    tuple val(patient), path(merged_vcf)

    output:
    tuple val(patient), path("${patient}")

    script:
    """
    mkdir ${patient}
    change_genome.pl ${patient} ${merged_vcf} ${params.bowtie_fa}
    bowtie2-build --threads ${task.cpus} ${patient}/${patient}.fa ${patient}
    """
}

process BOWTIE2_alt {
    cpus 15
    memory '20GB'
    maxForks 2

    publishDir "${params.outdir}/alt"

    input:
    tuple val(patient), val(sample), path(reads), path(alt_genome) 

    output:
    tuple val(sample), path("${sample}.alt.bam"), emit: bam
    tuple val(sample), path("${sample}.alt.bam.bai"), emit: bai

    script:
    def bowtie_cmd = reads.size() == 1 ? 
        "-U ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    bowtie2 --no-unal -p ${task.cpus} -x ${alt_genome}/${patient} ${bowtie_cmd} \\
        | samtools view -@ ${task.cpus} -b -q 10 -u - \\
        | samtools sort -@ ${task.cpus} -o ${sample}.alt.bam

    samtools index -@ ${task.cpus} ${sample}.alt.bam
    """
}

process DEMULTIPLEX_BAM {
    publishDir "${params.outdir}/sc_alignments"

    input:
    tuple val(sample), path(bam_file)

    output:
    path("*.bam"), emit: bam
    tuple path("*.bai"), emit: bai

    script:
    """
    for i in ${params.barcodes_dir}/${sample}*.txt; do
        subsample=\$(basename "\$i" .txt)

        if [[ ${bam_file} == *alt* ]]; then
            subsample=\${subsample}.alt
        else
            subsample=\${subsample}.ref
        fi

        samtools view -H ${bam_file} > SAM_header
        samtools view ${bam_file} | LC_ALL=C grep -F -f \$i > filtered_SAM_body
        cat SAM_header filtered_SAM_body | samtools view -b > \${subsample}.bam
        rm SAM_header filtered_SAM_body
    done

    samtools index -@ ${task.cpus} \${subsample}.bam
    """
}

process FILTER_VCF {
    maxForks 10
    cpus 1
    memory '3GB'

    publishDir "${params.outdir}/stats", mode: 'copy'

    input:
    tuple val(patient), val(sample), path(merged_bam), path(merged_vcf)

    output:
    path("*.stat")

    script:
    """
    filter_vcf.pl ${patient} ${sample} ${merged_bam} ${merged_vcf}
    """
}

process MERGE_VCF {
    publishDir "${params.outdir}/stats", mode: 'copy'

    input:
    path(stats)

    output:
    tuple path("ref.stat"), path("alt.stat")

    script:
    """
    merge_vcf.pl ${stats}
    """
}

process FINISH {
    publishDir "${params.outdir}/stats", mode: 'copy'

    input:
    path(stats)
    path(merged_vcfs)

    output:
    path("stat_nf.csv")

    script:
    """
    finish.pl ${stats} ${merged_vcfs}
    """
}

process CORRECT_BARCODES {
    publishDir "/media/leon/Polina/atac_rna/bc_corrected", mode: 'move'

    input:
    tuple path(fastq_file), val(assay_type)

    output:
    path("*.txt")

    script:
    """
    correct_barcode.R --barcode_file=${fastq_file} --assay_type=${assay_type}
    """
}

workflow {
    // Parsing the metadata table
    patient_dct = Channel.fromPath(params.metadata)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, row.patient) }
    
    // ATAC-seq
    ATAC_fastqs = Channel.of(50..83) // считывать sample_ids из таблицы с метадатой и искать в директории указанной в параметрах
        .map { n -> "SRR140487${n}" }
        .map { sample ->
            def fastq_f = "${params.atac_dir}/${sample}/${sample}_S1_L001_R1_001.fastq.gz"
            def fastq_b = "${params.atac_dir}/${sample}/${sample}_S1_L001_R2_001.fastq.gz"
            def fastq_r = "${params.atac_dir}/${sample}/${sample}_S1_L001_R3_001.fastq.gz"
        
            return [sample, [fastq_f, fastq_b, fastq_r], "ATAC"]
    }
    
    // Multiome ATAC-seq
    mmATAC_fastqs = Channel.of(351..398, 407..416)
        .map { n -> "SRR18593${n}" }
        .map { sample ->
            def fastq_f = "${params.mm_dir}/${sample}_1.fastq.gz"
            def fastq_b = "${params.mm_dir}/${sample}_2.fastq.gz"
            def fastq_r = "${params.mm_dir}/${sample}_3.fastq.gz"

            return [sample, [fastq_f, fastq_b, fastq_r], "mmATAC"]
        }

    // Multiome RNA-seq
    mmRNA_fastqs = Channel.of(339..350, 399..406)
        .map { n -> "SRR18593${n}" }
        .map { sample -> 
            def fastq_f = "${params.mm_dir}/${sample}_1.fastq.gz"
            def fastq_r = "${params.mm_dir}/${sample}_2.fastq.gz"

            return [sample, [fastq_f, fastq_r], "mmRNA"]
        }

    bc_correction_input = ATAC_fastqs.map { it -> it.flatten() }.map { it -> [it[2], it[4]] }
        .concat(mmATAC_fastqs.map { it -> it.flatten() }.map { it -> [it[2], it[4]] })
        .concat(mmRNA_fastqs.map { it -> it.flatten() }.map { it -> [it[1], it[3]] })
        
    CORRECT_BARCODES(bc_correction_input)
}

workflow rSNPs {
    // Parsing the metadata table
    patient_dct = Channel.fromPath(params.metadata)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, row.patient) }
    
    // ATAC-seq
    ATAC_fastqs = Channel.of(50..83) // считывать sample_ids из таблицы с метадатой и искать в директории указанной в параметрах
        .map { n -> "SRR140487${n}" }
        .map { sample ->
            def fastq_f = "${params.atac_dir}/${sample}/${sample}_S1_L001_R1_001.fastq.gz"
            def fastq_b = "${params.atac_dir}/${sample}/${sample}_S1_L001_R2_001.fastq.gz"
            def fastq_r = "${params.atac_dir}/${sample}/${sample}_S1_L001_R3_001.fastq.gz"
        
            return [sample, [fastq_f, fastq_b, fastq_r], "ATAC"]
    }
    
    // Multiome ATAC-seq
    mmATAC_fastqs = Channel.of(351..398, 407..416)
        .map { n -> "SRR18593${n}" }
        .map { sample ->
            def fastq_f = "${params.mm_dir}/${sample}_1.fastq.gz"
            def fastq_b = "${params.mm_dir}/${sample}_2.fastq.gz"
            def fastq_r = "${params.mm_dir}/${sample}_3.fastq.gz"

            return [sample, [fastq_f, fastq_b, fastq_r], "mmATAC"]
        }

    // Multiome RNA-seq
    mmRNA_fastqs = Channel.of(339..350, 399..406)
        .map { n -> "SRR18593${n}" }
        .map { sample -> 
            def fastq_f = "${params.mm_dir}/${sample}_1.fastq.gz"
            def fastq_r = "${params.mm_dir}/${sample}_2.fastq.gz"

            return [sample, [fastq_f, fastq_r], "mmRNA"]
        }

    bc_correction_input = ATAC_fastqs.map { it -> it.flatten() }.map { it -> [it[2], it[4]] }
        .concat(mmATAC_fastqs.map { it -> it.flatten() }.map { it -> [it[2], it[4]] })
        .concat(mmRNA_fastqs.map { it -> it.flatten() }.map { it -> [it[1], it[3]] })
        .view()
    CORRECT_BARCODES(bc_correction_input)

    // Preprocessing
    all_fastqs = ATAC_fastqs.concat(mmATAC_fastqs).concat(mmRNA_fastqs)
    trimmed_fastqs = PREPROCESS(all_fastqs).reads
    FASTQC(trimmed_fastqs)

    // Alignment to the reference genome
    ref_bams = BOWTIE2(trimmed_fastqs).bam

    // Merging all BAM files for each individual 
    merge_dct = ref_bams.join(patient_dct).map { it -> [it[2], it[1]] }.groupTuple() 
    merged_bams = MERGE_BAMS(merge_dct).bam

    // Identifying heterozygous positions 
    merged_vcfs = SNP_CALL(merged_bams)

    // Constructing an alternative genome
    alt_genomes = ALT_GENOME(merged_vcfs)

    // Alignment to the corresponding alternative genome
    bowtie_alt_input = trimmed_fastqs
        .join(patient_dct).map { it -> [it[2], it[0], it[1]] }
        .combine(alt_genomes, by: 0)
        
    alt_bams = BOWTIE_ALT(bowtie_alt_input).bam    
    
    // Demultipexing reference and alternative alignments by cell types
    sc_bams = DEMULTIPLEX_BAM(ref_bams.concat(alt_bams)).bam.flatten()
        .map { it -> [it.baseName, it] }

    // Calling SNPs in each BAM file separately and filtering against the merged BAM
    all_bams = ref_bams.concat(alt_bams).concat(sc_bams)

    filter_vcf_input = all_bams.map {it -> [it[0].split('_')[0], it[0], it[1]]}
        .combine(patient_dct, by: 0).map { it -> [it[3], it[0], it[1], it[2]] }
        .combine(merged_vcfs, by: 0).map {it -> [it[0], it[2], it[3], it[4]]}
    filtered_stats = FILTER_VCF(filter_vcf_input)

    // Collecting results
    merged_stats = MERGE_VCF(filtered_stats.collect())
    FINISH(merged_stats, merged_vcfs.map { it -> it[1] }. collect())
}

workflow asymmetry_shift {
   assay_dct = Channel.fromPath(params.metadata)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, row.assay_type) }

    fastqs = Channel.fromFilePairs(params.reads)
        .join(assay_dct)

    trimmed_fastqs = PREPROCESS(fastqs).reads
    FASTQC(trimmed_fastqs)

    alignments = HISAT2(trimmed_fastqs).bam

    GET_STATS(alignments)

    alignments.branch {
        rna: it[2] == 'RNASEQ'
        chip: it[2] == 'CHIPSEQ'
    }
    .set { assay }

    rna_alignments = assay.rna.map { sample, bam_path, assay_type -> bam_path }.collect()
    chip_alignments = assay.chip.map { sample, bam_path, assay_type -> bam_path }.collect()
    
    FEATURE_COUNTS(rna_alignments)
}

workflow {
    // rSNPs()
    asymmetry_shift()
}
