import sys
import re

from .logger import log

# Parsers for 'External_variation'

EXTERNAL = {
    'dbsnp_id': lambda x: x if x.startswith("rs") else None,
    'cosmic_id': lambda x: x if x.startswith("COSV") else None
}

def _get_external_ids(csq, EXTERNAL):
    ret = {k: "" for k in EXTERNAL}
    for i in csq['Existing_variation'].strip().split("&"):
        for key, func in EXTERNAL.items():
            val = func(i)
            if val:
                ret[key] += val if not ret[key] else ", " + val
    return ret

def get_csq_fields(vcf):
    desc = vcf.get_header_type("CSQ")['Description']
    ret = re.sub(".*Format: ", "", desc).split("|")
    return ret


def get_csq_dat(record, csq_fields, output_fields):
    ret = dict()
    for csq in record.INFO['CSQ'].split(","):
        dat = dict(zip(csq_fields, csq.split("|")))

        if not dat["ALLELE_NUM"]:
            log.error("Could not find 'ALLELE_NUM' in CSQ line")
            sys.exit(1)
        try:
            i = int(dat['ALLELE_NUM'])
        except:
            log.exception("Could not cast 'ALLELE_NUM' to int")
            sys.exit(1)

        if "PUBMED" in dat and dat["PUBMED"]:
            dat.update({"PUBMED": ", ".join(dat["PUBMED"].split("&"))})

        if "HGVSp" in dat:
            dat['HGVSp'] = dat['HGVSp'].replace("%3D", "=")

        ret[i] = {k: v for k, v in dat.items() if k in output_fields}
        ret[i].update(_get_external_ids(dat, EXTERNAL))
    return ret

