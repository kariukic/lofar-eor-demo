#! /usr/bin/env python3
import os
import casacore.tables as tab

def parse_ms_list(ms_list_or_txt):
    ms_list = [i for i in ms_list_or_txt if tab.tableexists(i)]
    if len(ms_list) == 0:
        ms_lists = []
        for ifile in ms_list_or_txt:
            if os.path.isfile(ifile):
                myf = open(ifile)
                ms_lists.append([i.strip() for i in myf if tab.tableexists(i.strip())])
    else:
        ms_lists = [ms_list]

    return [ms for ms_list in ms_lists for ms in ms_list]
