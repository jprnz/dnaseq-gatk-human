import sys
import yaml

from collections import OrderedDict
from .logger import log

def get_columns(args):
    try:
        columns = yaml.safe_load(open(args.config))['columns']
    except Exception:
        log.exception(f"Unable to load configuration file: {args.config}")
        sys.exit(1)

    return columns

def get_titles(columns):
    ret = list()
    for val in columns:
        keys = list(val.keys())
        if len(keys) != 1:
            log.error("Unexpected value in configuration file")
            sys.exit(1)
        ret += [val[keys[0]]['title']]
    return ret

def get_description(columns):
    ret = {'Column': list(), 'Description': list()}
    for val in [list(v.values())[0] for v in columns]:
            ret.update({
                'Column': ret['Column'] + [val['title']],
                'Description': ret['Description'] + [val['desc']]})
    return ret

