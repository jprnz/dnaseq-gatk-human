import os
import pandas as pd
import pickle

from .columns import get_description
from .logger import log
from .filter import filter
from collections import OrderedDict

EXCEL_LIMIT = 1048576 - 1

def _read_cache(dat, titles):
    chrom = list(dat.keys())[0]
    ret = pd.DataFrame(dat[chrom])
    ret.rename(columns=titles, inplace=True)
    ret = ret[titles.values()]
    return(chrom, ret)

def write_csv(csvfile, cache, columns):
    titles = OrderedDict([
        (list(v.keys())[0], v[list(v.keys())[0]]['title'])
        for v in columns
    ])

    first_chunk = True
    while True:
        if first_chunk:
            mode = 'w'
            header = True
        else:
            mode = 'a'
            header = False
        try:
            _, dat = _read_cache(pickle.load(cache), titles)
            dat.to_csv( csvfile, index=False, mode=mode, header=header)
            first_chunk = False
        except EOFError:
            break

def write_excel(dat, chrom, header, path):
    parts = list(range(0, len(dat), EXCEL_LIMIT)) + [len(dat)]

    ind_from = 0
    for part, ind_to in enumerate(parts[1:]):
        filename = f"{path}/{chrom}-{part:03}.xlsx"
        with pd.ExcelWriter(filename) as xlwrite:
            header.to_excel(xlwrite, index=False, sheet_name="Column Descriptions")
            dat.iloc[ind_from:ind_to].to_excel(xlwrite, index=False, sheet_name=f"{chrom}-{part:03}")
        ind_from = ind_to + 1

def write_excel_filtered(excelpath, cache, columns):
    header = pd.DataFrame(get_description(columns))
    titles = OrderedDict([
        (list(v.keys())[0], v[list(v.keys())[0]]['title'])
        for v in columns
    ])

    if not os.path.exists(excelpath):
        os.makedirs(excelpath)

    while True:
        try:
            chrom, dat = _read_cache(pickle.load(cache), titles)
            write_excel(filter(dat), chrom, header, excelpath)
        except EOFError:
            break
    

