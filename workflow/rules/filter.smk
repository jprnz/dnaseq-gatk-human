var_types = ["snp", "indel"]

# Determine output files
vqsr_targets = expand(filterdir + "/vqsr/{var}.vcf.gz", var=var_types)
hard_targets = expand(filterdir + "/hard/{var}.vcf.gz", var=var_types)

if filtering == "vqsr":
    var_targets = vqsr_targets
else:
    var_targets = hard_targets

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


rule merge_var_targets:
    input:
        var_targets
    output:
        vcf = filterdir + "/genotypes-filtered.vcf.gz",
        tbi = filterdir + "/genotypes-filtered.vcf.gz.tbi"
    log:
        filterdir + "/logs/merge-var-targets.log"
    conda:
        "../envs/gatk.yaml"
    params:
        vcfs = lambda wc, input: [f"-I {v}" for v in input]
    resources:
        mem_mb = 8000
    shell:
        "(set -x; "
        " gatk --java-options \"-Xms3G -Xmx8G -XX:ParallelGCThreads=5\""
        "   MergeVcfs"
        "   {params.vcfs}"
        "   --VALIDATION_STRINGENCY SILENT"
        "   --TMP_DIR {resources.tmpdir}"
        "   --OUTPUT {output.vcf}"
        ") &> {log}"

localrules: run_hard_filter
rule run_hard_filter:
    input: hard_targets

localrules: run_filter
rule run_filter:
    input: filterdir + "/genotypes-filtered.vcf.gz"


