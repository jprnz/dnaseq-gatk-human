rule bqsr_recal:
    input:
        bam = markdupdir + "/{sample}.bam",
        bai = markdupdir + "/{sample}.bai",
        ref_fasta = ref_fasta,
        ref_faidx = ref_fai,
        ref_dict = ref_dict,
        indels = known_indels,
        dbsnp = dbsnp_vcf
    output:
        bqsrdir + "/recal-files/{sample}.recal"
    log:
        bqsrdir + "/logs/recal-{sample}.log"
    params:
        regions = f"-L {regions}" if regions else "",
    conda:
        "../envs/gatk.yaml"
    resources:
        mem_mb = 10000
    shell:
        "(set -x; "
        " gatk --java-options \"-Xms8G -Xmx10G -XX:ParallelGCThreads=5 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10\""
        " BaseRecalibrator"
        "   --input {input.bam}"
        "   --reference {input.ref_fasta}"
        "   {params.regions}"
        "   --ip 500"
        "   --known-sites {input.dbsnp}"
        "   --known-sites {input.indels}"
        "   --tmp-dir {resources.tmpdir}"
        "   --output {output}"
        ") &> {log}"

rule bqsr_apply:
    input:
        bam = markdupdir + "/{sample}.bam",
        bai = markdupdir + "/{sample}.bai",
        recal = bqsrdir + "/recal-files/{sample}.recal",
        ref_fasta = ref_fasta
    output:
        bam = bqsrdir + "/{sample}.bam",
        bai = bqsrdir + "/{sample}.bai"
    log:
        bqsrdir + "/logs/apply-{sample}.log"
    params:
        regions = f"-L {regions}" if regions else "",
    conda:
        "../envs/gatk.yaml"
    resources:
        mem_mb = 16000
    shell:
        "(set -x; "
        " gatk --java-options \"-Xms8G -Xmx16G -XX:ParallelGCThreads=5 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10\""
        " ApplyBQSR"
        "   --input {input.bam}"
        "   --bqsr-recal-file {input.recal}"
        "   --reference {input.ref_fasta}"
        "   {params.regions}"
        "   --ip 500"
        "   --tmp-dir {resources.tmpdir}"
        "   --output {output.bam}"
        ") &> {log}"

localrules: run_bqsr
rule run_bqsr:
    input:
        expand(bqsrdir + "/{sample}.bam", sample=samples),
        expand(bqsrdir + "/recal-files/{sample}.recal", sample=samples)

