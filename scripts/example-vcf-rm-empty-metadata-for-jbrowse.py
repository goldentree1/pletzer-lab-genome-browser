#!/usr/bin/env python3

#
# USAGE:
# script.py < file.vcf > ready_for_jbrowse.vcf
#

import sys

for line in sys.stdin:
    sline = line.rstrip()
    if sline.startswith("##") and sline.find("=") == -1:
        pass
    else:
        print(sline)
