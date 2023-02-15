idxignore = "^GL*\|^KI*"

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
        "INPUT={input.bam} "
        "OUTPUT={output} "
        "REFERENCE_SEQUENCE={input.ref_fasta} "
        "INTERVALS={input.intervals} "
        "MINIMUM_MAPPING_QUALITY=1 "
        "MINIMUM_BASE_QUALITY=1 "
        "USE_FAST_ALGORITHM=false "
        "VALIDATION_STRINGENCY=SILENT) &> {log}"

#rule variant_metrics:
#    input:
#        intervals = regions
#        ref_fasta = ref_fasta
#    output:
#        metricsdir + "/wgs_metrics/{sample}.metrics"
#    log:
#        metricsdir + "/logs/wgs_metrics/{sample}.log"
#    conda:
#        "../envs/picard.yaml"
#    shell:
#        "(picard CollectVariantCallingMetrics "
#        "INPUT={input.bam} "
#        "OUTPUT={output} "
#        "REFERENCE_SEQUENCE={input.ref_fasta} "
#        "INTERVALS={input.intervals} "
#        "MINIMUM_MAPPING_QUALITY=1 "
#        "MINIMUM_BASE_QUALITY=1 "
#        "USE_FAST_ALGORITHM=false "
#        "VALIDATION_STRINGENCY=SILENT) &> {log}"

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

localrules: run_metrics
rule run_metrics:
    input:
        expand(metricsdir + "/wgs_metrics/{sample}.metrics", sample=samples),
        expand(metricsdir + "/idxstats/{sample}_idxstats.tsv", sample=samples),
        expand(metricsdir + "/samstats/{sample}_samstats.txt", sample=samples)

