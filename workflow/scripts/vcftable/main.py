import os
import sys
import tempfile

import pandas as pd

from cyvcf2 import VCF

from .logger import log
from .cli import get_args
from .csq import get_csq_fields
from .columns import get_columns
from .parser import parse_records
from .writer import write_csv, write_excel_filtered

try:
    TMPFILE = "test.pkl" #tempfile.NamedTemporaryFile()
except:
    log.exception("Could not create temporary file")
    sys.exit(1)

def vcftable():
    args = get_args()

    log.info(f"Reading config: {args.config}...")
    columns = get_columns(args)

    # Open VCF
    try:
        vcf = VCF(args.input)
    except:
        log.exception(f"Could not open VCF file at: {args.input}")
        sys.exit(1)

    # Set threads
    vcf.set_threads(args.threads)

    log.info(
        f"Prepairing to parse {len(columns)} fields of data "
        f"for {len(vcf.samples)} samples...")

    # Parse records
    log.info("Parsing CSV fields...")
    csq = get_csq_fields(vcf)
    with open(TMPFILE, 'wb+') as cache:
        parse_records(vcf, csq, columns, cache)

    # Write CSV file
    log.info("Writing CSV...")
    csvfile = f"{args.out_prefix}.csv"
    with open(TMPFILE, 'rb') as cache:
        write_csv(csvfile, cache, columns)

    # Write filtered per-chromosome / chunk
    log.info("Writing filtered excel file...")
    excelpath = f"{args.out_prefix}-VariantsByChromosome"
    with open(TMPFILE, 'rb') as cache:
        write_excel_filtered(excelpath, cache, columns)

