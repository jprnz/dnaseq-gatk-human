idxignore = "^GL*\|^KI*"

rule idxstats:
    input:
        bam =bqsrdir + "/{sample}.bam",
        bai =bqsrdir + "/{sample}.bai"
    output:
        metricsdir + "/idxstats/{sample}_idxstats.tsv"
    log:
        metricsdir + "/logs/idxstats/{sample}.log"
    params:
        ignore = idxignore
    conda:
        "../envs/samtools.yaml"
    shell:
        "(samtools idxstats {input.bam} | grep -v \"{params.ignore}\" > {output}) &> {log}"

rule samstats:
    input:
        bam = bqsrdir + "/{sample}.bam",
        bai = bqsrdir + "/{sample}.bai"
    output:
        metricsdir + "/samstats/{sample}_samstats.txt"
    log:
        metricsdir + "/logs/samstats/{sample}.log"
    params:
        ignore = idxignore
    conda:
        "../envs/samtools.yaml"
    shell:
        "(samtools stats {input.bam} > {output}) &> {log}"

rule wgs_metrics:
    input:
        bam = bqsrdir + "/{sample}.bam",
        bai = bqsrdir + "/{sample}.bai",
        intervals = regions,
        ref_fasta = ref_fasta
    output:
        metricsdir + "/wgs_metrics/{sample}.metrics"
    log:
        metricsdir + "/logs/wgs_metrics/{sample}.log"
    conda:
        "../envs/picard.yaml"
    shell:
        "(picard CollectWgsMetrics "
        "--INPUT {input.bam} "
        "--OUTPUT {output} "
        "--REFERENCE_SEQUENCE {input.ref_fasta} "
        "--INTERVALS {input.intervals} "
        "--MINIMUM_MAPPING_QUALITY 1 "
        "--MINIMUM_BASE_QUALITY 1 "
        "--USE_FAST_ALGORITHM false "
        "--VALIDATION_STRINGENCY SILENT) &> {log}"

rule variant_metrics:
    input:
        vcf = filterdir + "/genotypes-filtered.vcf.gz",
        intervals = regions,
        ref_fasta = ref_fasta,
        dbsnp = dbsnp_vcf,
    output:
        metricsdir + "/variant_metrics/filtered.variant_calling_detail_metrics",
        metricsdir + "/variant_metrics/filtered.variant_calling_summary_metrics"
    log:
        metricsdir + "/logs/variant_metrics/variant.log"
    conda:
        "../envs/picard.yaml"
    params:
        prefix = metricsdir + "/variant_metrics/filtered"
    threads: 10
    shell:
        "(picard CollectVariantCallingMetrics "
        "--INPUT {input.vcf} "
        "--DBSNP {input.dbsnp} "
        "--TARGET_INTERVALS {input.intervals} "
        "--OUTPUT {params.prefix} "
        "--REFERENCE_SEQUENCE {input.ref_fasta} "
        "--THREAD_COUNT {threads} "
        "--VALIDATION_STRINGENCY SILENT) &> {log}"

rule run_metrics:
    input:
        expand(metricsdir + "/idxstats/{sample}_idxstats.tsv", sample=samples),
        expand(metricsdir + "/samstats/{sample}_samstats.txt", sample=samples),
        expand(metricsdir + "/wgs_metrics/{sample}.metrics", sample=samples),
        metricsdir + "/variant_metrics/filtered.variant_calling_detail_metrics",
        metricsdir + "/variant_metrics/filtered.variant_calling_summary_metrics"

