#! /usr/bin/env python3

import os
import logging
import logging.config
import casacore.tables as tab
from argparse import ArgumentParser

from misc_utils import parse_ms_list

logging.config.fileConfig(os.path.join(os.path.dirname(__file__), 'logging.config'))
logging = logging.getLogger('root')


parser = ArgumentParser("get clip value from XY,YX rms of data")
parser.add_argument(
    "-i",
    "--ms_list",
    nargs="+",
    help="Measurement sets or file with list of MS",
    dest="ms_list",
    default="",
)


def main(args):
    logging.info("Checking measurement set(s)")
    ms_list = parse_ms_list(args.ms_list)
    logging.info(f"Total measurement sets found: {len(ms_list)}")

    for i in ms_list:
        logging.info("Flagging intrastrations")
        tab.taql(
            r'UPDATE $i SET FLAG=True WHERE mscal.baseline("/(.*)HBA0&\1HBA1/")'
        )  # mind the r


if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
