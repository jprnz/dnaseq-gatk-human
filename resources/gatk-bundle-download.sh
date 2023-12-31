#!/bin/bash
set -eo pipefail

url="https://storage.googleapis.com"
files=(
  "genomics-public-data/references/hg38/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf"
  "genomics-public-data/references/hg38/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf.idx"
  "genomics-public-data/references/hg38/v0/1000G_omni2.5.hg38.vcf.gz"
  "genomics-public-data/references/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi"
  "genomics-public-data/references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
  "genomics-public-data/references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi"
  "genomics-public-data/references/hg38/v0/1000G_phase3_v4_20130502.sites.hg38.vcf"
  "genomics-public-data/references/hg38/v0/1000G_phase3_v4_20130502.sites.hg38.vcf.idx"
  "genomics-public-data/references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz"
  "genomics-public-data/references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi"
  "genomics-public-data/references/hg38/v0/GRCh38.primary_assembly.genome.fa"
  "genomics-public-data/references/hg38/v0/GRCh38_gencode.v27.refFlat.txt"
  "genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.contam.UD"
  "genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.contam.V"
  "genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.contam.bed"
  "genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.contam.exome_calling_regions.v1.README.sh"
  "genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.contam.exome_calling_regions.v1.UD"
  "genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.contam.exome_calling_regions.v1.bed"
  "genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.contam.exome_calling_regions.v1.mu"
  "genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.contam.mu"
  "genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf"
  "genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx"
  "genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.dict"
  "genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta"
  "genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt"
  "genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.64.amb"
  "genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.64.ann"
  "genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.64.bwt"
  "genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.64.pac"
  "genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.64.sa"
  "genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.amb"
  "genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.ann"
  "genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.bwt"
  "genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
  "genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.pac"
  "genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.sa"
  "genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.haplotype_database.txt"
  "genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz"
  "genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"
  "genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.ref_cache.tar.gz"
  "genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.tile_db_header.vcf"
  "genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.vid"
  "genomics-public-data/references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
  "genomics-public-data/references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
  "genomics-public-data/references/hg38/v0/README"
  "genomics-public-data/references/hg38/v0/WholeGenomeShotgunContam.vcf"
  "genomics-public-data/references/hg38/v0/WholeGenomeShotgunContam.vcf.idx"
  "genomics-public-data/references/hg38/v0/autosomes-1kg-minusNA12878-ALL.vcf"
  #"genomics-public-data/references/hg38/v0/TileDB"
  #genomics-public-data/references/hg38/v0/CrossSpeciesContamination/ContaminantNormalizationFactors.txt"
  #genomics-public-data/references/hg38/v0/CrossSpeciesContamination/CrossSpeciesContaminant"
  #genomics-public-data/references/hg38/v0/CrossSpeciesContamination/CrossSpeciesContaminationMetrics.py"
  #genomics-public-data/references/hg38/v0/CrossSpeciesContamination/CrossSpeciesHuman"
)

for file in ${files[@]}; do
  wget -c "$url/$file"
done
