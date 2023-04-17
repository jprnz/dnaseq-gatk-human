import argparse

__doc__ = """\
Simple utility to convert a VCF to tabular csv/xlsx file.
"""

def get_args():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--input", "-i", metavar="vcf", required=True,
        help="VCF file to read")

    parser.add_argument(
        "--out-prefix", "-o", required=True,
        help="File extension will be added")

    parser.add_argument(
        "--xlsx", "-x", metavar="xlsx",
        help="Excel file to write")

    parser.add_argument(
        "--config", "-c", metavar="yaml",
        default="etc/config.yaml",
        help="Configuration file")

    parser.add_argument(
        "--threads", "-t", metavar="num", default=4,
        help="Number of threads for cyvcf2 to use")

    args = parser.parse_args()

    return args
