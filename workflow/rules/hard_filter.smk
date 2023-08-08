rule select_calls:
    input:
        haplotypedir + "/genotypes.vcf.gz"
    output:
        filterdir + "/hard/select-{var}.vcf.gz"
    log:
        filterdir + "/logs/hard-select-{var}.log"
    params:
        select = lambda wc: "SNP" if wc.var == "snp" else "INDEL"
    conda:
        "../envs/gatk.yaml"
    resources:
        mem_mb = 8000
    shell:
        "(set -x; "
        " gatk --java-options \"-Xms3G -Xmx8G -XX:ParallelGCThreads=5\""
        " SelectVariants"
        "   -V {input}"
        "   -select-type {params.select}"
        "   --tmp-dir {resources.tmpdir}"
        "   --output {output}"
        ") &> {log}"

rule hard_filter:
    input:
        filterdir + "/hard/select-{var}.vcf.gz"
    output:
        filterdir + "/hard/{var}.vcf.gz"
    log:
        filterdir + "/logs/hard-filter-{var}.log"
    conda:
        "../envs/gatk.yaml"
    params:
        filters = lambda wc: hard[wc.var]
    resources:
        mem_mb = 8000
    shell:
        "(set -x; "
        " gatk --java-options \"-Xms3G -Xmx8G -XX:ParallelGCThreads=5\""
        " VariantFiltration"
        "   -V {input}"
        "   --filter \"{params.filters}\""
        "   --filter-name HardFilter"
        "   --tmp-dir {resources.tmpdir}"
        "   --output {output}"
        ") &> {log}"

rule run_hard_filter:
    input: hard_targets

