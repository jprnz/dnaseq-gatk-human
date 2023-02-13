rule markdups:
    input:
        bam = bwadir + "/{sample}.bam",
        bai = bwadir + "/{sample}.bam.bai"
    output:
        bam = temp(markdupdir + "/{sample}.bam"),
        bai = temp(markdupdir + "/{sample}.bai"),
        metrics = markdupdir + "/metrics/{sample}.metrics"
    log:
        markdupdir + "/logs/{sample}.log"
    conda:
        "../envs/picard.yaml"
    resources:
        mem_mb = 8000
    shell:
        "(set -x; picard MarkDuplicates -Xmx8G "
        "  INPUT={input.bam} "
        "  OUTPUT={output.bam} "
        "  METRICS_FILE={output.metrics} "
        "  OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 "
        "  REMOVE_DUPLICATES=true "
        "  CREATE_INDEX=true "
        "  VALIDATION_STRINGENCY=LENIENT "
        ") &> {log}"

rule run_markdups:
    input:
        expand(markdupdir + "/{sample}.bam", sample=samples)

