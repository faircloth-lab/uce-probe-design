#!/usr/bin/env python
# encoding: utf-8
"""
File: rename_maf.py
Author: Brant Faircloth

Created by Brant Faircloth on 11 March 2012 19:03 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import glob
import argparse
from bx.align import maf

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Program description""")
    parser.add_argument(
            "indir",
            help="""The directory containing MAF files"""
        )
    parser.add_argument(
            "outdir",
            help="""The directory holding output files"""
        )
    parser.add_argument(
            "target",
            help="""The target prefix to add"""
        )
    parser.add_argument(
            "query",
            help="""The query prefix to add"""
        )
    return parser.parse_args()


def main():
    args = get_args()
    for f in glob.glob(os.path.join(args.indir, '*.maf')):
        inmaf = maf.Reader(open(f))
        outname = os.path.splitext(os.path.basename(f))[0] + ".rename.maf"
        outpth = os.path.join(args.outdir, outname)
        outf = maf.Writer(open(outpth, 'w'))
        for aln in inmaf:
            # change name
            aln.components[0].src = args.target + aln.components[0].src
            aln.components[1].src = args.query + aln.components[1].src
            outf.write(aln)
        outf.close()

if __name__ == '__main__':
    main()