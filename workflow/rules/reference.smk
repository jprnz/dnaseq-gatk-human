import os

from ftplib import FTP
from snakemake import shell
from operator import itemgetter


def get_ensembl_fasta(download=False, path=reference_path, log=None, **args):
    if not log:
        log = "/dev/null"

    datatype, species, build, release = itemgetter('datatype', 'species', 'build', 'release')(args)

    # Derived params
    branch = "grch37/" if release >= 81 and build == "GRCh37" else ""
    spec = f"{build}" if int(release) > 75 else f"{build}.{release}"
    species_cap = species.capitalize()

    # Define url
    url_base = f"ftp://ftp.ensembl.org/pub/{branch}release-{release}/fasta/{species}/{datatype}/{species_cap}.{spec}.{{suffix}}.gz"

    if datatype == "dna":
        # Check if primary_assembly is available
        suffix = "dna.primary_assembly.fa"
        url_test = url_base.format(suffix=suffix)
        try:
            shell("set -x; (curl --head -Ss -L {url_test}) &> {log}")
        except:
            suffix = "dna_sm.toplevel.fa"
    elif datatype == "cdna":
        suffix = "cdna.all.fa"
    elif datatype == "cds":
        suffix = "cds.all.fa"
    elif datatype == "ncrna":
        suffix = "ncrna.fa"
    elif datatype == "pep":
        suffix = "pep.all.fa"
    else:
        raise ValueError("invalid datatype, must be one of dna, cdna, cds, ncrna, pep")


    # Define paths
    local_path = f"{path}/{species}/{release}-{build}/{species_cap}.{spec}.{suffix}"
    url = url_base.format(suffix=suffix)

    if download:
        try:
            shell("set -x; (curl -S -L {url} | zcat > {local_path}) &> {log}")
        except:
            raise ValueError(
                    "Unable to download requested sequence data from Ensembl using urls:\n"
                    "{url}\nDid you check that this combination of species, build, "
                    "and release is actually provided?".format(url=url))
    else:
        return local_path


def get_ensembl_known(download=False, path=reference_path, log=None, **args):
    from snakemake import shell
    from operator import itemgetter


    species, build, release = itemgetter('species', 'build', 'release')(args)

    if species == "homo_sapiens":
        raise ValueError("Please use GATK bundle for human data")

    if not log:
        log = "/dev/null"

    # Get available file types
    ftp = FTP("ftp.ensembl.org")
    ftp.login()
    
    has_struc = any([
        True if "structural_variation" in v else False
        for v in ftp.nlst("/pub/release-{release}/variation/vcf/{species}/{species}")
        if v.endswith(".vcf.gz")])

    if has_struc:
        file_types = [".vcf.gz", "_structuarl_variation.vcf.gz"]
    else:
        file_types = [".vcf.gz"]

    if download:
        base_url = f"ftp://ftp.ensembl.org/pub/release-{release}/variation/vcf/{species}/{species}"

        for file_type in file_types:
            url = base_url + file_type
            local = f"{path}/{species}/{release}-{build}/{species}{file_type}"
            try:
                shell("set -x; (curl -L {url} > {local}) &> {log}")
            except Exception as e:
                raise ValueError(
                    "Unable to download requested sequence data from Ensembl using url:\n"
                    "{url}\n Did you check that this combination of species, build, "
                    "and release is actually provided?".format(url=url))
    else:
        ret = [f"{path}/{species}/{release}-{build}/{species}{v}" for v in file_types]
        return ret


genome_fasta = get_ensembl_fasta(datatype="dna", **reference_args)
known_snps = get_ensembl_known(**reference_args)
known_snps_idx = [v + ".tbi" for v in known_snps]
genome_dict = os.path.splitext(genome_fasta)[0] + ".dict"
genome_faidx = genome_fasta + ".fai"

rule download_genome:
    output:
        genome_fasta
    log:
        reference_path + "/logs/{species}-{release}-{build}-dna-fasta.log".format(**reference_args)
    run:
        get_ensembl_fasta(datatype='dna', download=True, log=log[0], **reference_args)

rule download_known:
    output:
        known_snps
    log:
        reference_path + "/logs/{species}-{release}-{build}-dna-known.log".format(**reference_args)
    run:
        get_ensembl_known(download=True, log=log[0], **reference_args)

rule known_index:
    input:
        known_snps
    output:
        known_snps_idx
    log:
        reference_path + "/logs/{species}-{release}-{build}-dna-known-index.log".format(**reference_args)
    conda:
        "../envs/gatk.yaml"
    shell:
        '(set -x; for fn in "{input}"; do gatk IndexFeatureFile -I $fn; done) &> {log}'

rule genome_faidx:
    input:
        genome_fasta
    output:
        f"{genome_fasta}.fai"
    log:
        reference_path + "/logs/{species}-{release}-{build}-dna-fai.log".format(**reference_args)
    conda:
        "../envs/samtools.yaml"
    shell:
        "(set -x; samtools faidx {input}) &> {log}"

rule genome_dict:
    input:
        genome_fasta
    output:
        genome_dict
    log:
        reference_path + "/logs/{species}-{release}-{build}-dna-dict.log".format(**reference_args)
    conda:
        "../envs/picard.yaml"
    shell:
        "(set -x; picard CreateSequenceDictionary "
        "  --REFERENCE {input}"
        "  --OUTPUT {output}) &> {log}"

rule bwa_index:
    input:
        genome_fasta
    output:
        bwa_index + ".sa"
    log:
        reference_path + "/logs/{species}-{release}-{build}-dna-bwa.log".format(**reference_args)
    params:
        prefix = bwa_index
    conda:
        "../envs/bwa.yaml"
    shell:
        "(set -x; bwa index "
        "  -p {params.prefix} "
        "  {input}) &> {log}"

#rule star_index_genome:
#    input:
#        fasta = genome_fasta,
#        gtf = genome_gtf
#    output:
#        directory(star_index)
#    log:
#        star_index + ".log"
#    conda:
#        "../envs/star.yaml"
#    resources:
#        mem_mb = 48000
#    threads: 16
#    shell:
#        "(set -x; STAR"
#        "  --genomeDir {output} "
#        "  --runMode genomeGenerate "
#        "  --genomeFastaFiles {input.fasta} "
#        "  --sjdbGTFfile {input.gtf} "
#        "  --sjdbOverhang 100 "
#        "  --runThreadN {threads} "
#        ") &> {log}"


