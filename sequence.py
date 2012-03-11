#!/usr/bin/env python
# encoding: utf-8
"""
sequence.py

Created by Brant Faircloth on 2008-07-08.
Copyright (c) 2008 Brant Faircloth. All rights reserved.

This module is meant to be used with the summary.py program for MAF files.
"""
import pdb
import numpy
import re

class stack(numpy.ndarray):
    '''Takes a sequence stack (ndarray of strings) and performs several operations'''
    def __new__(subtype, data, info='Sequence Stack', dtype='|S1', copy=False):
        assert data.ndim > 1, 'sequence is not a stack'
        # Make sure we are working with an array, and copy the data if requested
        subarr = numpy.array(data, dtype=dtype, copy=copy)
        # Transform 'subarr' from an ndarray to our new subclass.
        subarr = subarr.view(subtype)
        # Use the specified 'info' parameter if given
        if info is not None:
            subarr.info = info
        # Otherwise, use data's info attribute if it exists
        elif hasattr(data, 'info'):
            subarr.info = data.info
        # Finally, we must return the newly created object:
        return subarr
    
    def __array_finalize__(self,obj):
        # We use the getattr method to set a default if 'obj' doesn't have the 'info' attribute
        self.info = getattr(obj, 'info', {})
    
    def __repr__(self):
        desc = """array(data=%(data)s,tag='%(tag)s')"""
        return desc % {'data': str(self), 'tag':self.info }
    
    def iupac(self, s):
        '''define translation for mismatched bases, forgetting gaps'''
        iupacCode = {
        'A':'A',
        'C':'C',
        'G':'G',
        'T':'T',
        'AG':'R',
        'CT':'Y',
        'AC':'M',
        'GT':'K',
        'AT':'W',
        'CG':'S',
        'CGT':'B',
        'AGT':'D',
        'ACT':'H',
        'ACG':'V',
        'ACGT':'N'
        }
        return iupacCode[s]
    
    def consensus(self):
        '''loops over consensus => dumb consensus and maintaining mask status'''
        # ensure that sequence is indeed a stack
        assert self.ndim > 1, 'sequence is not a stack'
        # setup regex to search for lowercase bases (masking)
        lc = re.compile('[acgt]+')
        cons = []
        # for i in the length of the row
        for b in range(len(self[0,:])):
            # slice the data across columns (i,e. each colum is a
            # stack of bases at the aligned position).
            stringRep = self[:,b]
            # if gap anywhere, leave it and only it
            if '-' in stringRep:
                cons.append('-')
            # there appear to be some 'N'|'n' bases.  If in the sequence column,
            # represent all bases as N|n, maintaining case.
            elif 'N' in stringRep:
                cons.append('N')
            elif 'n' in stringRep:
                cons.append('n')
            # otherwise, see if there are mismatches and translate
            else:
                bases = []
                # ignore same base > 1 time, uppercase all
                [bases.append(i.upper()) for i in stringRep if i.upper() not in bases]
                # sort the bases, so they match the dict entries
                bases.sort()
                # translate, if needed
                try:
                    base = self.iupac(''.join(bases))
                except KeyError:
                    pdb.set_trace()
                # if the column contains a masked base, from any 
                # species, lowercase the base
                string = ''.join(stringRep)
                if lc.search(string):base = base.lower()
                cons.append(base)
        # return the consensus as a string
        self.cons = ''.join(cons)
        return self.cons
            