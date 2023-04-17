import re
import sys

from .logger import log

def get_variant_dat(record):
    ret = dict()
    pos = {
            "chrom": record.CHROM,
            "pos": record.POS,
            "ref": record.REF}
    for i, alt in enumerate(record.ALT, start=1):
        ret[i] = {
            "chrom": record.CHROM,
            "pos": record.POS,
            "ref": record.REF,
            "alt": alt}
    ret.update({0: ret[1].copy()})
    ret[0]['alt'] = None
    return ret







