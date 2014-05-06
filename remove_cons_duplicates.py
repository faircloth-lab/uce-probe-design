#!/usr/bin/env python
# encoding: utf-8

"""
remove_fish_cons_duplicates.py

Created by Brant Faircloth on 29 April 2010 11:54 PDT (-0700).
Copyright (c) 2010 Brant C. Faircloth. All rights reserved.

PURPOSE:  This file takes an output file from lastz and locates sequences for
which there are duplicate matches (e.g. matches to something other than self).
The script then updates the `cons` database in mysql by adding a binary/boolean
column and filling it for those conserved regions that are duplicates.

USAGE:  python remove_fish_cons_duplicates.py --configuration=db.conf \
            --input=my_lastz_file.lastz
"""

import pdb
import os
import sys
import sqlite3
import optparse



def interface():
    '''Command-line interface'''
    usage = "usage: %prog [options]"

    p = optparse.OptionParser(usage)

    p.add_option('--db', action='store', 
type='string', default = None, help='The path to the cons database.', 
metavar='FILE')
    p.add_option('--lastz', dest = 'input', action='store', 
type='string', default = None, help='The path to the configuration file.', 
metavar='FILE')

    (options,arg) = p.parse_args()
    return options, arg 


def main():
    options, arg = interface()
    conn = sqlite3.connect(options.db)
    cur = conn.cursor()
    lastz = open(options.input, 'rU')
    all_ids = set()
    #dupes = {}
    dupe_set = set()
    for line in lastz:

        line_elements = line.strip('\n').split('\t')
        name1, name2 = line_elements[1].strip('>'), line_elements[6].strip('>')
        id1, id2 = name1.split('_')[0], name2.split('_')[0]
        all_ids.add(id1)
        all_ids.add(id2)
        if id1 == id2:
            pass
        else:
            dupe_set.add(id1)
            dupe_set.add(id2)
    #pdb.set_trace()
    non_dupes = all_ids.difference(dupe_set)
    try:
        cur.execute('ALTER TABLE cons ADD COLUMN duplicate int(1)')
    except sqlite3.OperationalError, e:
        if e == "duplicate column name: duplicate":
            cur.execute('UPDATE cons set duplicate = NULL')
    for d in dupe_set:
        cur.execute('UPDATE cons set duplicate = 1 WHERE id = ?', (d,))
    for nd in non_dupes:
        cur.execute('UPDATE cons set duplicate = 0 WHERE id = ?', (nd,))
    lastz.close()
    cur.close()
    conn.commit()
    conn.close()



if __name__ == '__main__':
    main()