#!/usr/bin/env python

import os
import sys
import json
import glob

EXE = "python ../../server/problem2txt.py"
BONUS_TYPES = ["", "GLOBALIST", "BREAK_A_LEG", "WALLHACK"]


def convert_bonus_format(infile, bonus_type, outfile):
    os.system("{} {} {} > {}".format(EXE, infile, bonus_type, outfile))


def main(args):
    if len(args) < 2:
        print("python3 {} infiles...".format(args[0]))
        sys.exit()

    infiles = []
    for arg in args[1:]:
        infiles.extend(glob.glob(arg))

    for infile in infiles:
        for bonus_type in BONUS_TYPES:
            print("infile: {}, bonus: {}".format(infile, bonus_type))
            outfile = infile.replace("problem", "in")
            if bonus_type:
                outfile = "{}.{}".format(outfile, bonus_type)
            convert_bonus_format(infile, bonus_type, outfile)


if __name__ == "__main__":
    main(sys.argv)
