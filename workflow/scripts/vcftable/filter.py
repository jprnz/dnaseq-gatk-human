import sys

import pandas as pd

from .logger import log

def filter(dat):
    return dat[
        (
            dat['gnomAD AF Exome'].notna() 
            | dat['gnomAD AF Genome'].notna() 
            | dat['dbSNP'].notna()
        ) & (pd.to_numeric(dat['Total read count']) > 5)
    ]

