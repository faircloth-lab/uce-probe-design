#!/usr/bin/env python
# encoding: utf-8
"""
File: fasta_renamer.py
Author: Brant Faircloth

Created by Brant Faircloth on 12 March 2012 15:03 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import argparse
import subprocess
import pdb

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Helps mask a contig, preserving fasta header""")
    parser.add_argument(
            "infile",
            help="""The input file to mask"""
        )
    parser.add_argument(
            "--name-portion",
            dest='name',
            type=int,
            default=-1,
            help="""The portion of the fasta header to use as the new name"""
        )
    parser.add_argument(
            "--prefix",
            dest='prefix',
            type=str,
            default="",
            help="""A prefix to add to the new chromo names"""
        )
    parser.add_argument(
            "--no-mask",
            dest='no_mask',
            action="store_true",
            default=False,
            help="""The portion of the fasta header to use as the new name"""
        )
    return parser.parse_args()


def main():
    #pdb.set_trace()
    args = get_args()
    names = {}
    temp1 = "{}.temp1".format(args.infile)
    temp2 = "{}.temp2".format(args.infile)
    outf = open(temp1, 'w')
    mask_file = os.path.splitext(args.infile)[0] + ".rm.out"
    if args.no_mask:
        for line in open(args.infile, 'rU'):
            if line.startswith('>'):
                print line
                newname = args.prefix + line.lstrip('>').split(',')[0].split(' ')[args.name]
                newline = ">{}\n".format(newname)
                outf.write(newline)
            else:
                outf.write(line)
    else:
        for line in open(args.infile, 'rU'):
            if line.startswith('>'):
                print line
                gb = line.split('|')[3]
                newname = args.prefix + line.split(',')[0].split(' ')[args.name]
                names[gb] = newname
                newline = ">{}\n".format(gb)
                outf.write(newline)
            else:
                outf.write(line)
        outf.close()
        cmd = ["maskOutFa", "-softAdd", temp1, mask_file, temp2]
        subprocess.Popen(cmd).wait()
        final = "{}.masked".format(args.infile)
        outf = open(final, 'w')
        for line in open(temp2, 'rU'):
            #pdb.set_trace()
            if line.startswith('>'):
                #pdb.set_trace()
                iden = line.strip('>').strip()
                line = ">{}\n".format(names[iden])
                print line
                outf.write(line)
            else:
                outf.write(line)
        try:
            os.remove(temp1)
            os.remove(temp2)
        except:
            pass
    outf.close()




if __name__ == '__main__':
    main()
