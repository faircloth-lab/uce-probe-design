#!/usr/bin/env python
# encoding: utf-8
"""
summary.py

Created by Brant Faircloth on 2008-07-05.
Copyright (c) 2008 Brant Faircloth. All rights reserved.

This program scans MAF files for conserved elements and stores 
those results in an sqlite database

"""

import os
import re
import sys
import pdb # remove at some point
import time
import numpy
import MySQLdb
import sequence
import optparse
import ConfigParser
import bx.align.maf
import multiprocessing


def interface():
    '''Get the starting parameters from a configuration file'''
    usage = "usage: %prog [options]"

    p = optparse.OptionParser(usage)

    p.add_option('--configuration', dest = 'conf', action='store', \
type='string', default = None, help='The path to the configuration file.', \
metavar='FILE')
    p.add_option('--maf', dest = 'maf', action='store', \
type='string', default = None, help='The path to the directory containing maf file(s).', \
metavar='FILE')
    p.add_option('--alignment-length', dest = 'align', action='store', \
type='int', default = 25, help='The minimum acceptable alignment length.')
    p.add_option('--consensus-length', dest = 'consensus', action='store', \
type='int', default = 25, help='The minimum acceptable consensus length.')
    p.add_option('--metadata-key', dest = 'metadata', action='store', \
type='string', default = 25, help='The _primary_ species in the alignment \
(e.g. the one on top).')

    (options,arg) = p.parse_args()
    if not options.conf or not options.metadata:
        p.print_help()
        sys.exit(2)
    if not os.path.isfile(options.conf):
        print "You must provide a valid path to the configuration file."
        p.print_help()
        sys.exit(2)
    return options, arg

def spScreen(a, minAlignLength):
    '''screen alignments to ensure minSpecies and minAlignLength'''
    for spp in a.components:
        if len(a.components[0].text) > minAlignLength:
            return a

def alignMetadata(counter, candAlign, cons, refPosition, altPosition, metadataKey):
    '''get metdata for alignment based on species in metadataKey'''
    for seq in candAlign.components:
        name = seq.src
        metadata = {}
        #pdb.set_trace()
        if name.split('.')[0] == metadataKey:
            metadata['target_spp']      = name.split('.')[0]
            metadata['target_chromo']   = name.split('.')[1]
            metadata['target_start']    = seq.forward_strand_start
            metadata['target_len']      = seq.size
            metadata['target_end']      = seq.forward_strand_end
            metadata['target_strand']   = seq.strand
            metadata['cons']            = cons
            metadata['cons_len']        = len(cons)
            # add values to metadata, making up for 0 indexing
            metadata['target_cons_start']      = metadata['target_start'] + 1 + refPosition[0]
            metadata['target_cons_end']        = metadata['target_start'] + refPosition[1]
            metadata['query_spp']       = candAlign.components[1].src.split('.')[0]
            metadata['query_chromo']    = candAlign.components[1].src.split('.')[1]
            metadata['query_strand']    = candAlign.components[1].strand
            metadata['query_len']       = candAlign.components[1].size
            # deal with forward and reverse strand weirdness
            #pdb.set_trace()
            if metadata['query_strand'] == '+':
                metadata['query_start'] = candAlign.components[1].start
                metadata['query_end'] = candAlign.components[1].start + candAlign.components[1].size
                metadata['query_cons_start'] = candAlign.components[1].start + 1 + altPosition[0]
                metadata['query_cons_end']   = candAlign.components[1].start + altPosition[1]
            else:
                metadata['query_end'] = candAlign.components[1].src_size - candAlign.components[1].start
                metadata['query_start'] = metadata['query_end'] - (candAlign.components[1].size - 1)
                metadata['query_cons_end']   = candAlign.components[1].src_size - (candAlign.components[1].start + altPosition[0])
                metadata['query_cons_start'] = candAlign.components[1].src_size - (candAlign.components[1].start + altPosition[1] - 1)
            metadata['target_cons_map']      = (('%s:%s-%s') % (metadata['target_chromo'], metadata['target_cons_start'], metadata['target_cons_end']))
            metadata['query_cons_map']       = (('%s:%s-%s') % (metadata['query_chromo'], metadata['query_cons_start'], metadata['query_cons_end']))
            #metadata['seq']             = (('%s_%s') % (metadata['target_chromo'], counter))
        break
    #pdb.set_trace()
    return metadata

def createCons(candAlign):
    '''stack sequence and return dumb (but smart!) consensus with
    metadata'''
    for seq in range(len(candAlign.components)):
        if seq == 0:
            zString = candAlign.components[seq].text
            zString = numpy.array(list(zString))
            seqArray = zString
        else:
            nzString = candAlign.components[seq].text
            nzString = numpy.array(list(nzString))
            seqArray = numpy.vstack((seqArray,nzString))
    #pdb.set_trace()
    seqStack = sequence.stack(seqArray)
    consensus = seqStack.consensus()
    return consensus

def filterCons(unfilteredConsensus, minConsensusLength, iterate = False):
    '''filter out alignments with short, gappy, mismatching shit (most of them)'''
    # find masked|unmasked block > minConsensusLength
    searchString = (('[ACGT]{%i,}') % (minConsensusLength))
    pattern = re.compile(searchString)
    if not iterate:
        masked = pattern.search(unfilteredConsensus)
        if masked:
            return list(masked.group())
        else:
            return False
    else:
        masked = pattern.findall(unfilteredConsensus)
        if masked:
            return masked
        else:
            return False
        #return masked.group()
    #else:
    #    return False

def positioner(candAlign, cons):
    '''return correct positions of the conserved area relative to the reference seq
    by degapping while also dealing with repeat-masked sequence in the conserved area'''
    # strip gap character from reference seq
    pattern = re.compile('-+')
    cleanCandAlign = pattern.sub('', candAlign.text)
    # deal with upper/lowercase issues btw reference <--> alt and
    # repeat-masked bases
    caseUnawareCons = []
    for letter in cons:
        if letter.isupper():
            bracket = (('[%s%s]') % (letter, letter.lower()))
            caseUnawareCons.append(bracket)
        else:
            bracket = (('[%s%s]') % (letter, letter.upper()))
            caseUnawareCons.append(bracket)
    caseUnawareCons = ''.join(caseUnawareCons)
    # find position of conserved sequence relative to gapless
    # candAlign
    pattern = re.compile(caseUnawareCons)
    position = pattern.search(cleanCandAlign)
    return position.span()
    
def createConsTable(cur):
    '''create a table to hold the results'''
    try:
        # if previous tables exist, drop them
        # TODO: fix createDbase() to drop tables safely
        cur.execute('''DROP TABLE cons''')
    except:
        pass
    # create the primers results table
    cur.execute('''CREATE TABLE cons (
    id MEDIUMINT UNSIGNED NOT NULL AUTO_INCREMENT,
    target_spp VARCHAR(7) NOT NULL,
    target_chromo VARCHAR(20) NOT NULL,
    target_start INT UNSIGNED NOT NULL,
    target_end INT UNSIGNED NOT NULL,
    target_len SMALLINT UNSIGNED NOT NULL,
    target_strand varchar(1) NOT NULL,
    target_cons_start INT UNSIGNED NOT NULL,
    target_cons_end INT UNSIGNED NOT NULL,
    target_cons_map VARCHAR(30) NOT NULL,
    query_spp VARCHAR(7) NOT NULL,
    query_chromo VARCHAR(20) NOT NULL,
    query_start INT UNSIGNED NOT NULL,
    query_end INT UNSIGNED NOT NULL,
    query_len SMALLINT UNSIGNED NOT NULL,
    query_cons_start INT UNSIGNED NOT NULL,
    query_cons_end INT UNSIGNED NOT NULL,
    query_cons_map VARCHAR(30) NOT NULL,
    query_strand varchar(1) NOT NULL,
    cons TEXT NOT NULL,
    cons_len SMALLINT UNSIGNED NOT NULL,
    duplicate BOOLEAN,
    PRIMARY KEY (id),
    INDEX cons_cons_len(cons_len),
    INDEX cons_target_chromo(target_chromo),
    INDEX cons_query_chromo(query_chromo)
    ) ENGINE=InnoDB
    ''')
    
def store(cur, metadata):
    '''store the results in mysql'''
    cur.execute('''insert into cons (
    target_spp, 
    target_chromo, 
    target_start, 
    target_end,
    target_len, 
    target_strand,
    target_cons_start,
    target_cons_end, 
    target_cons_map,
    query_spp,
    query_chromo,
    query_start,
    query_end,
    query_len,
    query_cons_start,
    query_cons_end,
    query_cons_map,
    query_strand,
    cons,
    cons_len) 
    values 
    (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)''', 
    (metadata['target_spp'],
    metadata['target_chromo'],
    metadata['target_start'],
    metadata['target_end'],
    metadata['target_len'],
    metadata['target_strand'],
    metadata['target_cons_start'],
    metadata['target_cons_end'],
    metadata['target_cons_map'],
    metadata['query_spp'],
    metadata['query_chromo'],
    metadata['query_start'],
    metadata['query_end'],
    metadata['query_len'],
    metadata['query_cons_start'],
    metadata['query_cons_end'],
    metadata['query_cons_map'],
    metadata['query_strand'],
    metadata['cons'], 
    metadata['cons_len']))

def worker(input, minConsensusLength, minAlignLength, metadataKey, conf):
    # we need a separate connection for each mysql cursor or they are going
    # start going into locking hell and things will go poorly. Creating a new 
    # connection for each worker process is the easiest/laziest solution.
    # Connection pooling (DB-API) didn't work so hot, but probably because 
    # I'm slightly retarded.
    conn = MySQLdb.connect(
    	user=conf.get('Database','USER'), 
    	passwd=conf.get('Database','PASSWORD'), 
    	db=conf.get('Database','DATABASE')
		)
    cur = conn.cursor()
    file = open(input,'rU')
    parser = bx.align.maf.Reader(file)
    a = parser.next()
    # select only those alignments of > minSpecies
    counter = 0
    META = {}
    while a:
        counter += 1
        candAlign = spScreen(a, minAlignLength)
        if candAlign:
            # create sequence stack and stack -> dumb consensus
            unfilteredConsensus = createCons(candAlign)
            # filter out consensi with < 1 contiguous block of minConsensus
            conserved = filterCons(unfilteredConsensus, minConsensusLength, True)
            #pdb.set_trace()
            if conserved:
                for cons in conserved:
                    #print '%s: ****Valid consensus****' % counter
                    #print cons
                    # find 'real' positions in reference sequence (galGal3 here)
                    # by degapping
                    refPosition = positioner(candAlign.components[0], cons)
                    # find 'real' positions in alternate sequence (anoCar1 here)
                    # by degapping
                    altPosition = positioner(candAlign.components[1], cons)
                    # get sequence metadata
                    metadata = alignMetadata(counter, candAlign, cons, refPosition, altPosition, metadataKey)
                    # store start, totalLength, end, consensus somewhere
                    # insert records to dbase
                    store(cur, metadata)
                    conn.commit()  
        a = parser.next()
    # close the MAF reader
    parser.close()
    # close the file
    file.close()
    # keep our connection load low
    cur.close()
    conn.close()

def file_gen(directory):
    '''create an iterable list of filenames in the appropriate directory'''
    for f in os.listdir(directory):
        if os.path.splitext(f)[1] == '.maf' and os.path.isfile(os.path.join(directory, f)):
            yield os.path.join(directory, f)

def main():
    start = time.time()
    conf = ConfigParser.ConfigParser()
    options, arg = interface()
    conf.read(options.conf)
	#minConsensusLength = 40
    #minAlignLength     = 40
    #metadataKey        = 'taeGut1'
    #directory = '/Users/bcf/data/genome/taeGut1/maf'
    n_procs = multiprocessing.cpu_count() - 2
    #n_procs = 1
    # connect to our dbase
    conn = MySQLdb.connect(
        user=conf.get('Database','USER'), 
        passwd=conf.get('Database','PASSWORD'), 
        db=conf.get('Database','DATABASE')
    )
    #minConsensusLength = 25
    #minAlignLength     = 25
    #metadataKey        = 'oryLat2'
    #directory = '/Users/bcf/Data/alignments/oryLat2.danRer6/maf'
    #n_procs = multiprocessing.cpu_count() - 2
    n_procs = 1
    # connect to our dbase
    cur = conn.cursor()
    # create a new table or drop if exists
    createConsTable(cur)
    conn.commit()
    # get an iterator of *.maf files from directory
    files = file_gen(options.maf)
    #pdb.set_trace()
    if n_procs > 1:
        print 'Multiprocessing.  Number of processors = ', n_procs
        try:
            threads = []
            while files:
                if len(threads) < n_procs:
                    p = multiprocessing.Process(target=worker, args=(files.next(), options.consensus, options.align, options.metadata, conf))
                    p.start()
                    threads.append(p)
                else:
                    for p in threads:
                        if not p.is_alive():
                            threads.remove(p)
        except StopIteration:
            pass
    
    else:
        print 'Not using multiprocessing'
        try:
            while files:
                worker(files.next(), options.consensus, options.align, options.metadata, conf)
        except StopIteration:
            pass
    # commit any remaining changes
    conn.commit()
    cur.close()
    conn.close()
    # finish up execution time
    end = time.time()
    execution = (end - start)/60.
    print 'Time for execution = %f min.' % (execution)
    
if __name__ == '__main__':
    main()