import sys
import yaml
import pandas as pd

from .logger import log

EXCEL_LIMIT = 1048576 - 5

def _write_excel_file(dat, header, filename):
    with pd.ExcelWriter(filename) as xlwriter:
        header.to_excel(xlwriter, sheet_name="Column Descriptions", index=False)
        dat.to_excel(xlwriter, sheet_name="Variants", float_format="%.4f", index=False)

def write_excel(dat, header, filename):
    if len(dat) > EXCEL_LIMIT:
        parts = list(range(EXCEL_LIMIT, len(dat), EXCEL_LIMIT))
        for part in parts:
            import ipdb; ipdb.set_trace()
    else:
        _write_excel_file(dat, header, filename)
    

