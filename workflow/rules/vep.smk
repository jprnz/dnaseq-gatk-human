vep_url = f"https://ftp.ensembl.org/pub/release-{vep_version}/variation/indexed_vep_cache/homo_sapiens_vep_{vep_version}_GRCh38.tar.gz"
vep_cache = vep_cache_path + f"/homo_sapiens/{vep_version}_GRCh38"
vep_cache_target = vep_cache + "/info.txt"


rule vep_cache:
  output: 
    vep_cache_target
  log:
    vep_cache_path + f"/logs/homo_sapiens-GRCh38-{vep_version}.log"
  params:
    path = vep_cache_path,
    filename = os.path.basename(vep_url)
  shell:
    "(set -x;"
    " cd {params.path} "
    " && wget -nq {vep_url} "
    " && tar vxf {params.filename}"
    ") &> {log}"

rule vep:
  input:
    vcf = filterdir + "/genotypes-filtered.vcf.gz",
    ref = ref_fasta,
    cache = vep_cache_target,
  output:
    vcf = vepdir + "/genotypes-filtered-annotated.vcf.gz",
    htm = vepdir + "/vep-summary.html"
  log:
    vepdir + "/logs/vep.log"
  params:
    cache = vep_cache_path
  conda:
    "../envs/vep.yaml"
  resources:
    mem_mb = 16000
  threads: 10
  shell:
    "(set -x; vep "
    "  --input_file {input.vcf}"
    "  --dir_cache {params.cache}"
    "  --fasta {input.ref}"
    "  --format vcf"
    "  --verbose"
    "  --cache"
    "  --allele_number"
    "  --everything"
    "  --pick_allele"
    "  --force"
    "  --vcf"
    "  --buffer_size 50000"
    "  --compress_output gzip"
    "  --output_file {output.vcf}"
    "  --stats_file {output.htm}"
    "  --fork {threads}"
    ") &> {log}"

rule run_vep:
  input: rules.vep.output
