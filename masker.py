#!/usr/bin/env python
# encoding: utf-8
"""
File: masker.py
Author: Brant Faircloth

Created by Brant Faircloth on 11 March 2012 15:03 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import argparse
import subprocess
from seqtools.sequence import fasta

import pdb

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Helps mask a contig, preserving fasta header""")
    parser.add_argument(
            "infile",
            help="""The input file to mask"""
        )
    return parser.parse_args()


def main():
    #pdb.set_trace()
    args = get_args()
    names = {}
    temp1 = "{}.temp1".format(args.infile)
    temp2 = "{}.temp2".format(args.infile)
    outf = fasta.FastaWriter(temp1)
    mask_file = os.path.splitext(args.infile)[0] + ".fa.out"
    f = fasta.FastaReader(args.infile)
    for seq in f:
        print seq.identifier
        gb = seq.identifier.split('|')[3]
        newname = seq.identifier.split(',')[0].split(' ')[-1]
        names[gb] = newname
        seq.identifier = ">{}".format(gb)
        outf.write(seq)
    outf.close()
    cmd = ["maskOutFa", "-softAdd", temp1, mask_file, temp2]
    subprocess.Popen(cmd).wait()
    final = "{}.masked".format(args.infile)
    outf = fasta.FastaWriter(final)
    for seq in fasta.FastaReader(temp2):
        iden = seq.identifier.strip('>')
        seq.identifier = "{}".format(names[iden])
        print seq.identifier
        outf.write(seq)
    outf.close()




if __name__ == '__main__':
    main()
