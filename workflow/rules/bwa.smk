rule bwa:
    input:
        r1 = fastpdir + "/{sample}_R1.fastq.gz",
        r2 = fastpdir + "/{sample}_R2.fastq.gz",
        index = ref_bwt
    output:
        bam = temp(bwadir + "/{sample}.bam"),
        bai = temp(bwadir + "/{sample}.bam.bai")
    log:
        bwadir + "/logs/{sample}.log"
    conda:
        "../envs/bwa.yaml"
    params:
        extra = r"-R '@RG\tID:{sample}\tSM:{sample}\tPU:{sample}\tLB:Illumina\tPL:Illumina'",
        index = bwa_index
    resources:
        mem_mb = 16000
    threads: 10
    shell:
        "(set -x; "
        "  bwa mem "
        "    -t {threads} "
        "    {params.extra} "
        "    {params.index} "
        "    {input.r1} "
        "    {input.r2} "
        "  | samtools sort -o {output.bam} "
        "  && samtools index -b {output.bam} "
        ") &> {log}"

localrules: run_bwa
rule run_bwa:
    input:
        expand(bwadir + "/{sample}.bam", sample=samples)

