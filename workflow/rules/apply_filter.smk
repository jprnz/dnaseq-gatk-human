# Wildcards 
var_types = ["snp", "indel"]

# Targets
hard_targets = expand(filterdir + "/hard/{var}.vcf.gz", var=var_types)
vqsr_targets = expand(filterdir + "/vqsr/{var}.vcf.gz", var=var_types)

# Source rules and get targets
if filtering == "vqsr":
    var_targets = vqsr_targets
    include: "vqsr_filter.smk"
else:
    var_targets = hard_targets
    include: "hard_filter.smk"

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

rule run_filter:
    input: rules.merge_var_targets.output

