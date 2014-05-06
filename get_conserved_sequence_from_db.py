#!/usr/bin/env python
# encoding: utf-8
"""
getConservedFishSequence.py

Created by Brant Faircloth on 2010-04-28.
Copyright (c) 2010 Brant C. Faircloth. All rights reserved.

PURPOSE: This program filters sequence (from ./summary.py) stored in a mysql and
returns a fasta file of those sequences making it past the filter.

USAGE: python getConservedFishSequence.py --configuration=db.conf
"""


import pdb
import sqlite3
import argparse
import bx.seq.twobit


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Program description""")
    parser.add_argument(
            "db",
            help="""The cons database"""
        ),
    parser.add_argument(
            "output",
            help="""The output file"""
        )
    parser.add_argument(
            "--conserved-length",
            dest="cons",
            type=int,
            default=40,
            help="""The minimum conserved length to return""",
        ),
    parser.add_argument(
            "--buffer-to",
            dest="buffer",
            type=int,
            default=None,
            help="""The flanking sequence to return""",
        ),
    parser.add_argument(
            "--two-bit",
            dest="twobit",
            type=str,
            default=None,
            help="""The path to the primary (target) twobit file.""",
        ),
    parser.add_argument(
            "--drop-dupes",
            dest="dupes",
            action="store_true",
            default=False,
            help="""Drop duplicate sequences""",
        )
    return parser.parse_args()


def sequence_writer(handle, iden, species, chromo, start, end, sequence):
    header = '{0}_{1}_{2}_{3}_{4}'.format(iden, species, chromo, start, end)
    handle.write('>{0}\n{1}\n'.format(header, sequence))


def get_reads_drop_dupes(cur, cons):
    cur.execute("""
                SELECT id,
                target_spp,
                target_chromo,
                target_cons_start,
                target_cons_end,
                cons
                FROM cons
                WHERE duplicate != 1
                AND cons_len >= ?""", (cons,)
            )
    return cur.fetchall()


def get_reads(cur, cons):
    cur.execute("""
                SELECT id,
                target_spp,
                target_chromo,
                target_cons_start,
                target_cons_end,
                cons
                FROM cons
                WHERE cons_len >= ?""", (cons,)
            )
    return cur.fetchall()


def main():
    args = get_args()
    #pdb.set_trace()
    conn = sqlite3.connect(args.db)
    cur = conn.cursor()
    if args.dupes:
        data = get_reads_drop_dupes(cur, args.cons)
    else:
        data = get_reads(cur, args.cons)
    handle = open(args.output, 'w')
    if args.buffer:
        tb = bx.seq.twobit.TwoBitFile(file(args.twobit))
    for d in data:
        iden, target_spp, target_chromo, target_cons_start, \
            target_cons_end, cons = d
        if args.buffer and len(cons) < args.buffer:
            # determine the difference in length
            diff = abs(len(cons) - args.buffer)
            # split the difference to each side of the conserved
            if diff % 2:
                diff = (diff + 1) / 2
            else:
                diff /= 2
            start = target_cons_start - diff
            end = target_cons_end + diff
            if start < 0:
                start, end = 0, args.buffer
            # get buffered sequence
            #pdb.set_trace()
            cons_buff = tb[target_chromo][start:end]
            sequence_writer(handle, iden, target_spp, target_chromo, start, \
                end, cons_buff)
        else:
            sequence_writer(handle, iden, target_spp, target_chromo, \
                target_cons_start, target_cons_end, cons)
    #pdb.set_trace()
    handle.close()
    cur.close()
    conn.close()

if __name__ == '__main__':
    main()
