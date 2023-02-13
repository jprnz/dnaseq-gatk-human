rule bwa:
    input:
        r1 = fastpdir + "/{sample}_R1.fastq.gz",
        r2 = fastpdir + "/{sample}_R2.fastq.gz",
        index = bwa_index + ".sa"
    output:
        bam = temp(bwadir + "/{sample}.bam"),
        bai = temp(bwadir + "/{sample}.bam.bai")
    log:
        bwadir + "/logs/{sample}.log"
    conda:
        "../envs/bwa.yaml"
    params:
        extra = r"-R '@RG\tID:{sample}\tSM:{sample}\tPU:{sample}\tLB:Illumina\tPL:Illumina'",
        index = bwa_index.split(".sa")[0],
    resources:
        mem_mb = 8000
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

rule run_bwa:
    input:
        expand(bwadir + "/{sample}.bam", sample=samples)

