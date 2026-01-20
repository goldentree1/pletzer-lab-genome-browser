#!/usr/bin/env python3

import sys

"""
Replaces GFF lines with Name=PLES... with the old_locus_tag= attribute, because
apparently that's how they should be named.

Use like a bash cmd (stdin/stdout).

TESTER:
    cat test/LESB58_ASM2664v1.gff.sorted.gff.noregion.renamedPLES | wc -l
    cat test/LESB58_ASM2664v1.gff.sorted.gff.noregion | wc -l

^ firstly, make sure those output same n of lines.

TESTER2:
    cat scripts/gff_ples_name_corrector_TEST1.txt | python scripts/GFF_PLES_name_corrector.py | diff scripts/gff_ples_name_corrector_TEST1.txt -
    # should be 3 diff'd thingys above, cuz there were 3 Name=PLES files.
"""

for line in sys.stdin:
    sline = line.rstrip()
    fields = sline.split("\t")
    if len(fields) == 9:
        attributes = fields[8]
        if "Name=PLES" in attributes and "old_locus_tag=" in attributes:
            nameAttr = ""
            oldLocusTagAttr = ""
            attrList = []
            for a in attributes.split(";"):
                if a.startswith("Name=PLES"):
                    nameAttr = a.replace("Name=", "")
                elif a.startswith("old_locus_tag="):
                    oldLocusTagAttr = a.replace("old_locus_tag=", "")
                    attrList.append(a)
                else:
                    attrList.append(a)
            attrList.insert(1, f"Name={oldLocusTagAttr}")
            newAttributesString = ";".join(attrList)
            print(
                f"{fields[0]}\t{fields[1]}\t{fields[2]}\t{fields[3]}\t{fields[4]}\t{fields[5]}\t{fields[6]}\t{fields[7]}\t{newAttributesString}"
            )
        else:
            print(sline)
    else:
        print(sline)
