if scatter_split_intervals:
    mode = "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW"
else:
    mode = "INTERVAL_SUBDIVISION"

if haplotype_bamfiles:
    haplotype_bam = "--bam-output " + haplotypedir + "/scatter/{sample}/{scatter}.bam"
else:
    haplotype_bam = False

checkpoint rule make_intervals:
    input:
        regions = regions,
        ref_fasta = ref_fasta
    output:
        directory(haplotypedir + "/intervals")
    log:
        haplotypedir + "/logs/make-intervals.log"
    params:
        scatter_count = scatter_count,
        mode = mode
    conda:
        "../envs/gatk.yaml"
    resources:
        mem_mb = 3000
    shell:
        "(set -x; "
        "  gatk --java-options \"-Xms3000m -Xmx3250m\""
        "  SplitIntervals"
        "    --input {input.regions}"
        "    --reference {input.ref_fasta}"
        "    --scatter-count {params.scatter_count}"
        "    --subdivision-mode {params.mode}"
        "    --interval-merging-rule OVERLAPPING_ONLY"
        "    --tmp-dir {resources.tmpdir}"
        "    --output {output}"
        ") &> {log}"

def haplotype_gvcfs(wildcards):
    sample = wildcards.get("sample", samples)
    happath = haplotypedir + "/scatter/{sample}/{scatter}.g.vcf.gz"
    invpath = checkpoints.make_intervals.get(**wildcards).output[0]
    scatter = glob_wildcards(invpath + "/{i}-scattered.interval_list")
    ret = expand(happath, sample=sample, scatter=scatter.i)
    return sorted(ret)

def haplotype_bams(wildcards):
    sample = wildcards.get("sample", samples)
    happath = haplotypedir + "/scatter/{sample}/{scatter}.bam"
    invpath = checkpoints.make_intervals.get(**wildcards).output[0]
    scatter = glob_wildcards(invpath + "/{i}-scattered.interval_list")
    ret = expand(happath, sample=sample, scatter=scatter.i)
    return sorted(ret)

def get_scatter_targets(wildcards):
    ret = [haplotypedir + "/gvcf_files/{sample}.g.vcf.gz"]
    if haplotype_bamfiles:
        ret.append(haplotypedir + "/bam_files/{sample}.bam")
    return ret

rule haplotype_caller:
    input:
        bam = bqsrdir + "/{sample}.bam",
        bai = bqsrdir + "/{sample}.bai",
        interval = haplotypedir + "/intervals/{scatter}-scattered.interval_list",
        ref_fasta = ref_fasta,
        dbsnp = dbsnp_vcf
    output:
        vcf = haplotypedir + "/scatter/{sample}/{scatter}.g.vcf.gz",
        bam = haplotypedir + "/scatter/{sample}/{scatter}.bam",
        bai = haplotypedir + "/scatter/{sample}/{scatter}.bai"
    log:
        haplotypedir + "/logs/scatter/{sample}/{scatter}.log"
    params:
        bam_output = haplotype_bam if haplotype_bam else ""
    conda:
        "../envs/gatk.yaml"
    resources:
        mem_mb = 12000
    shell:
        "(set -x; "
        " gatk --java-options \"-Xms8G -Xmx10G -XX:ParallelGCThreads=5 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10\""
        " HaplotypeCaller"
        "   --input {input.bam}"
        "   --reference {input.ref_fasta}"
        "   --intervals {input.interval}"
        "   --dbsnp {input.dbsnp}"
        "   -G StandardAnnotation"
        "   -G AS_StandardAnnotation"
        "   -ERC GVCF"
        "   {params.bam_output}"
        "   --tmp-dir {resources.tmpdir}"
        "   --output {output.vcf}"
        " && touch {output.bam} {output.bai}) &> {log}"

rule gather_haplotypes:
    input:
        haplotype_gvcfs
    output:
        vcf = haplotypedir + "/gvcf_files/{sample}.g.vcf.gz",
    log:
        haplotypedir + "/logs/{sample}-merge-vcfs.log"
    conda:
        "../envs/gatk.yaml"
    params:
        vcfs = lambda wc, input: [f"-I {v}" for v in input]
    resources:
        mem_mb = 4000
    shell:
        "(set -x; "
        " gatk --java-options \"-Xms2G -Xmx3G -XX:ParallelGCThreads=5\""
        " MergeVcfs"
        "   {params.vcfs}"
        "   --VALIDATION_STRINGENCY SILENT"
        "   --TMP_DIR {resources.tmpdir}"
        "   --OUTPUT {output}"
        ") &> {log}"

rule gather_bamfiles:
    input:
        haplotype_bams
    output:
        bam = haplotypedir + "/bam_files/{sample}.bam"
    log:
        haplotypedir + "/logs/{sample}-merge-bams.log"
    conda:
        "../envs/gatk.yaml"
    params:
        bams = lambda wc, input: [f"-I {v}" for v in input]
    resources:
        mem_mb = 4000
    shell:
        "(set -x; "
        " gatk --java-options \"-Xms2G -Xmx3G -XX:ParallelGCThreads=5\""
        " MergeSamFiles"
        "   {params.bams}"
        "   --CREATE_INDEX"
        "   --VALIDATION_STRINGENCY SILENT"
        "   --TMP_DIR {resources.tmpdir}"
        "   --OUTPUT {output.bam}"
        ") &> {log}"
        
rule remove_scatter:
    input:
        get_scatter_targets
    output:
        haplotypedir + "/logs/{sample}-remove-scatter.log"
    params:
        path = haplotypedir + "/scatter/{sample}"
    priority: 10
    shell:
        "(set -x; rm -rv {params.path}) &> {output}"

rule combine_gvcfs:
    input:
        vcfs = expand(haplotypedir + "/gvcf_files/{sample}.g.vcf.gz", sample=samples),
        intervals = regions
    output:
        directory(haplotypedir + "/genotype_db")
    log:
        haplotypedir + "/logs/genotypedb.log"
    conda:
        "../envs/gatk.yaml"
    params:
        vcfs = lambda wc, input: [f"-V {v}" for v in input.vcfs]
    resources:
        mem_mb = 25000
    threads: 5
    shell:
        "(set -x; "
        " gatk --java-options \"-Xms8000m -Xmx25000m -XX:ParallelGCThreads=5\""
        " GenomicsDBImport"
        "   {params.vcfs}"
        "   --genomicsdb-workspace-path {output}"
        "   --intervals {input.intervals}"
        "   --merge-input-intervals"
        "   --reader-threads 5"
        "   --tmp-dir {resources.tmpdir}"
        ") &> {log}"

rule genotype_gvcfs:
    input:
        gendb = haplotypedir + "/genotype_db",
        intervals = regions,
        ref_fasta = ref_fasta,
        dbsnp = dbsnp_vcf
    output:
        haplotypedir + "/genotypes.vcf.gz"
    log:
        haplotypedir + "/logs/genotype.log"
    conda:
        "../envs/gatk.yaml"
    resources:
        mem_mb = 25000
    threads: 1
    shell:
        "(set -x; "
        " gatk --java-options \"-Xms8000m -Xmx25000m -XX:ParallelGCThreads=5\""
        " GenotypeGVCFs"
        "   -V gendb://{input.gendb}"
        "   --reference {input.ref_fasta}"
        "   --intervals {input.intervals}"
        "   --dbsnp {input.dbsnp}"
        "   -G StandardAnnotation"
        "   -G AS_StandardAnnotation"
        "   --merge-input-intervals"
        "   --tmp-dir {resources.tmpdir}"
        "   --output {output}"
        ") &> {log}"


hap_targets = [haplotypedir + "/genotypes.vcf.gz"]
bam_targets = expand(haplotypedir + "/bam_files/{sample}.bam", sample=samples)

if haplotype_bamfiles:
    hap_targets += bam_targets  

localrules: remove_scatter, run_haplotype
rule run_haplotype:
    input: 
        expand(haplotypedir + "/logs/{sample}-remove-scatter.log", sample=samples),
        haplotypedir + "/genotypes.vcf.gz"

