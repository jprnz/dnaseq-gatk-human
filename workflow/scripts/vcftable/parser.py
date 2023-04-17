import pickle

from rich.progress import track

from .logger import log
from .csq import get_csq_dat
from .variant import get_variant_dat
from .sample import get_sample_dat

EXCEL_LIMIT = 1048576 - 1

def _get_allele_feilds(i, allele, var, csq, fields):
    ret = allele

    if i in var:
        ret.update(**var[i])
    else:
        # allele i is -1 (not called)
        pass

    if i in csq:
        ret.update(**csq[i])
    else:
        # allele i is 0 or -1
        pass
    ret.update({k: None for k in fields if k not in ret}) 
    return ret


def _parse_record(record, samples, csq_fields, out_fields):
    csq = get_csq_dat(record, csq_fields, out_fields)
    var = get_variant_dat(record)
    ind = get_sample_dat(record, samples)

    # Iterate over samples, then alleles
    ret = list()
    for sample in ind.values():
        for i, allele in sample.items():
            if i == 0 and allele['gt'] != '0/0':
                continue
            dat = _get_allele_feilds(i, allele, var, csq, out_fields)
            ret.append(dat)

    return ret


def parse_records(vcf, csq_fields, columns, cache):

    # Parse columns
    out_fields = sum([list(k.keys()) for k in columns], [])

    # Samples
    samples = vcf.samples

    # Keep track of current chromosome
    this_chrom = None

    ret = list()
    for record in track(vcf, description="Processing VCF: "):
        if not this_chrom:
            this_chrom = record.CHROM

        dat = _parse_record(record, samples, csq_fields, out_fields)
        if (this_chrom != record.CHROM):
            pickle.dump({this_chrom: ret}, cache)
            this_chrom = record.CHROM
            ret = dat
        else:
            ret.extend(dat)

    # In case of single chromosome
    pickle.dump({record.CHROM: ret}, cache)
    return ret


