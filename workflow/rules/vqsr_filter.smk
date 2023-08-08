rule sites_only_vcf:
    input:
        haplotypedir + "/genotypes.vcf.gz"
    output:
        filterdir + "/vqsr/sites-only.vcf.gz"
    log:
        filterdir + "/logs/sites-only.log"
    conda:
        "../envs/gatk.yaml"
    resources:
        mem_mb = 8000
    shell:
        "(set -x; "
        " gatk --java-options \"-Xms3G -Xmx8G -XX:ParallelGCThreads=5\""
        " MakeSitesOnlyVcf"
        "   -I {input}"
        "   -O {output}"
        "   --TMP_DIR {resources.tmpdir}"
        ") &> {log}"

rule vqsr_snp_recal:
    input:
        vcf = filterdir + "/vqsr/sites-only.vcf.gz",
        hapmap = hapmap_vcf,
        omni = omni_vcf,
        onekg = onekg_vcf,
        dbsnp = dbsnp_vcf
    output:
        recal = filterdir + "/vqsr/snp.recal",
        tranches = filterdir + "/vqsr/snp.tranches"
    log:
        filterdir + "/logs/vqsr-snp-recal.log"
    params:
        ann = [f"-an {v} " for v in vqsr['snp']['annotations']],
        tranches = [f"-tranche {v} " for v in vqsr['snp']['tranches']]
    conda:
        "../envs/gatk.yaml"
    resources:
        mem_mb = 8000
    shell:
        "(set -x; "
        " gatk --java-options \"-Xms3G -Xmx8G -XX:ParallelGCThreads=5\""
        " VariantRecalibrator"
        "   -V {input.vcf}"
        "   -mode SNP"
        "   --use-allele-specific-annotations"
        "   --max-gaussians 6"
        "   {params.ann}"
        "   {params.tranches}"
        "   --trust-all-polymorphic"
        "   --resource:hapmap,known=false,training=true,truth=true,prior=15 {input.hapmap}"
        "   --resource:omni,known=false,training=true,truth=true,prior=12 {input.omni}"
        "   --resource:1000G,known=false,training=true,truth=false,prior=10 {input.onekg}"
        "   --resource:dbsnp,known=true,training=false,truth=false,prior=7 {input.dbsnp}"
        "   --tranches-file {output.tranches}"
        "   --tmp-dir {resources.tmpdir}"
        "   --output {output.recal}"
        ") &> {log}"

rule vqsr_indel_recal:
    input:
        vcf = filterdir + "/vqsr/sites-only.vcf.gz",
        hapmap = hapmap_vcf,
        mills = mills_vcf,
        axiom = axiom_vcf,
        dbsnp = dbsnp_vcf
    output:
        recal = filterdir + "/vqsr/indel.recal",
        tranches = filterdir + "/vqsr/indel.tranches"
    log:
        filterdir + "/logs/vqsr-indel-recal.log"
    params:
        ann = [f"-an {v} " for v in vqsr['indel']['annotations']],
        tranches = [f"-tranche {v} " for v in vqsr['indel']['tranches']]
    conda:
        "../envs/gatk.yaml"
    resources:
        mem_mb = 8000
    shell:
        "(set -x; "
        " gatk --java-options \"-Xms3G -Xmx8G -XX:ParallelGCThreads=5\""
        " VariantRecalibrator"
        "   -V {input.vcf}"
        "   -mode INDEL"
        "   --use-allele-specific-annotations"
        "   --max-gaussians 3"
        "   {params.ann}"
        "   {params.tranches}"
        "   --trust-all-polymorphic"
        "   -resource:mills,known=false,training=true,truth=true,prior=12 {input.mills}"
        "   -resource:axiomPoly,known=false,training=true,truth=false,prior=10 {input.axiom}"
        "   -resource:dbsnp,known=true,training=false,truth=false,prior=2 {input.dbsnp}"
        "   --tranches-file {output.tranches}"
        "   --tmp-dir {resources.tmpdir}"
        "   --output {output.recal}"
        ") &> {log}"

rule vqsr_apply:
    input:
        vcf = haplotypedir + "/genotypes.vcf.gz",
        snp_tranches = filterdir + "/vqsr/snp.tranches",
        snp_recal = filterdir + "/vqsr/snp.recal",
        indel_tranches = filterdir + "/vqsr/indel.tranches",
        indel_recal = filterdir + "/vqsr/indel.recal"
    output:
        vcf = filterdir + "/genotypes-filtered.vcf.gz",
        idx = filterdir + "/genotypes-filtered.vcf.gz.tbi",
        tmp = temp(filterdir + "/vqsr/snp.tmp.vcf.gz"),
        tmp_idx = temp(filterdir + "/vqsr/snp.tmp.vcf.gz.tbi"),
    log:
        filterdir + "/logs/vqsr-apply.log"
    conda:
        "../envs/gatk.yaml"
    params:
        snp_level = vqsr['snp']['filter_level'],
        indel_level = vqsr['indel']['filter_level'],
    resources:
        mem_mb = 8000
    shell:
        "(set -x; "
        " gatk --java-options \"-Xms3G -Xmx8G -XX:ParallelGCThreads=5\""
        " ApplyVQSR"
        "   -V {input.vcf}"
        "   -mode SNP"
        "   --recal-file {input.snp_recal}"
        "   --tranches-file {input.snp_tranches}"
        "   --truth-sensitivity-filter-level {params.snp_level}"
        "   --use-allele-specific-annotations"
        "   --tmp-dir {resources.tmpdir}"
        "   --output {output.tmp} &&"
        " gatk --java-options \"-Xms3G -Xmx8G -XX:ParallelGCThreads=5\""
        " ApplyVQSR"
        "   -V {output.tmp}"
        "   -mode INDEL"
        "   --recal-file {input.indel_recal}"
        "   --tranches-file {input.indel_tranches}"
        "   --truth-sensitivity-filter-level {params.indel_level}"
        "   --use-allele-specific-annotations"
        "   --tmp-dir {resources.tmpdir}"
        "   --output {output.vcf}"
        ") &> {log}"

rule run_vqsr:
    input: rules.vqsr_apply.output

