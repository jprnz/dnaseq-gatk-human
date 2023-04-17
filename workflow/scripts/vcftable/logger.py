#!/bin/env python
import logging

from rich.logging import RichHandler

logging.basicConfig(
    level="NOTSET",
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(
        rich_tracebacks=True,
        show_path=False,
        omit_repeated_times=False)])

log = logging.getLogger("rich")



