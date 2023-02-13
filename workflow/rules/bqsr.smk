rule bqsr_recal:
    input:
        bam = markdupdir + "/{sample}.bam",
        bai = markdupdir + "/{sample}.bai",
        genome_fasta = genome_fasta,
        genome_faidx = genome_faidx,
        genome_dict = genome_dict,
        known = known_snps,
        known_idx = known_snps_idx
    output:
        bqsrdir + "/recal-files/{sample}.recal"
    log:
        bqsrdir + "/logs/recal-{sample}.log"
    params:
        regions = "-L f{regions_file}" if regions else "",
        known = [f"--known-sites {v}" for v in known_snps]
    conda:
        "../envs/gatk.yaml"
    resources:
        mem_mb = 8000
    shell:
        "(set -x; "
        "  gatk BaseRecalibrator"
        "   --input {input.bam}"
        "   --reference {input.genome_fasta}"
        "   {params.regions}"
        "   {params.known}"
        "   --output {output}"
        ") &> {log}"

rule bqsr_apply:
    input:
        bam = markdupdir + "/{sample}.bam",
        bai = markdupdir + "/{sample}.bai",
        recal = bqsrdir + "/recal-files/{sample}.recal",
        genome_fasta = genome_fasta
    output:
        bam = bqsrdir + "/{sample}.bam",
        bai = bqsrdir + "/{sample}.bai"
    log:
        bqsrdir + "/logs/apply-{sample}.log"
    params:
        regions = "-L f{regions_file}" if regions else ""
    conda:
        "../envs/gatk.yaml"
    resources:
        mem_mb = 8000
    shell:
        "(set -x; "
        "  gatk ApplyBQSR"
        "   --input {input.bam}"
        "   --bqsr-recal-file {input.recal} "
        "   --reference {input.genome_fasta} "
        "   {params.regions} "
        "   --output {output.bam} "
        ") &> {log}"

rule run_bqsr:
    input:
        expand(bqsrdir + "/{sample}.bam", sample=samples),
        expand(bqsrdir + "/recal-files/{sample}.recal", sample=samples)

