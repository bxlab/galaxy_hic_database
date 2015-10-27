#!/usr/bin/env python

import sys
import os
import h5py
import argparse

import numpy

def main():
    parser = generate_parser()
    args = parser.parse_args()
    Manager( args.fend, args.in_prefix )
    return None

def generate_parser():
    """Generate an argument parser."""
    description = "%(prog)s -- Create a raw file of paired aligned reads for a HiC experiment from fastq or sam files."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument( dest="fend", type=str, action='store', help="HiFive fend file name." )
    parser.add_argument( dest="in_prefix", type=str, action='store', help="File prefix for paired reads." )
    return parser

class Manager( object ):

    def __init__( self, fend_fname, in_prefix ):
        self.load_fends( fend_fname )
        self.find_num_partitions( in_prefix )
        for i in range(self.num_partitions):
            pairer = Pairer( self.fends, "%s_1.sam.%i" % (in_prefix, i),
                             "%s_2.sam.%i" % (in_prefix, i), self.seq_len )
            pairer.pair_reads()

    def load_fends( self, fname ):
        self.fends = {}
        infile = h5py.File( fname, 'r' )
        chr_indices = infile[ 'chr_indices' ][ ... ]
        chromosomes = infile[ 'chromosomes' ][ ... ]
        starts = infile[ 'fends' ][ 'start' ][ ... ]
        stops = infile[ 'fends' ][ 'stop' ][ ... ]
        for i in range( chr_indices.shape[ 0 ] - 1 ):
            self.fends[ chromosomes[ i ] ] = numpy.r_[ starts[ chr_indices[ i ]:chr_indices[ i + 1 ]:2 ],
                                                       stops[ chr_indices[ i + 1 ] - 1 ] ]
        infile.close()
        return None    

    def find_num_partitions( self, in_prefix ):
        for line in open("%s_1.stats" % in_prefix):
            temp = line.rstrip('\n').split(': ')
            if temp[0] == 'Total reads':
                self.num_partitions = (int(temp[1].strip(' ')) - 1) / 2000000 + 1
            elif temp[0] == 'Sequence length':
                self.seq_len = int(temp[1].strip(' '))

class Pairer( object ):

    def __init__( self, fends, in_fname1, in_fname2, seq_len ):
        self.fends = fends
        self.mapper1 = Reader( in_fname1 )
        self.mapper2 = Reader( in_fname2 )
        self.raw_fname = in_fname1.replace('_1.sam', '.raw')
        self.output = open(self.raw_fname, 'w')
        self.stats_fname = in_fname1.replace('_1.sam', '.stats')
        self.paired = 0
        self.invalid_chimeric = 0
        self.valid_chimeric = 0
        self.valid = 0
        self.seq_len = seq_len
        return None

    def print_stats( self ):
        output = open(self.stats_fname, 'w')
        print >> output, "Paired reads: %i" % (self.paired)
        print >> output, "Valid paired reads: %i" % (self.valid)
        print >> output, "Total chimeric reads: %i" % (self.invalid_chimeric + self.valid_chimeric)
        print >> output, "Valid chimeric reads: %i" % (self.valid_chimeric)
        output.close()

    def pair_reads( self ):
        for key, read1 in self.mapper1.reads.iteritems():
            if key in self.mapper2.reads:
                self.paired += 1
                read2 = self.mapper2.reads[key]
                if read1 is None or read2 is None:
                    self.invalid_chimeric += 1
                elif len(read1) > 4 or len(read2) > 4:
                    self.resolve_chimeric( read1, read2 )
                else:
                    self.valid += 1
                    self.resolve_standard( read1[:3], read2[:3] )
        self.print_stats()

    def resolve_standard( self, read1, read2 ):
        if read1[2] == '-':
            read1[1] = str(int(read1[1]) + self.seq_len)
        if read2[2] == '-':
            read2[1] = str(int(read2[1]) + self.seq_len)
        if read1 < read2:
            print >> self.output, '\t'.join(read1 + read2)
        else:
            print >> self.output, '\t'.join(read2 + read1)

    def resolve_chimeric( self, read1, read2 ):
        valid = True
        if len(read1) > 4:
            fend1, read1 = self.find_paired_fends( read1 )
            if len(read2) > 4:
                fend2, read2 = self.find_paired_fends( read2 )
                if fend1 != -self.find_fend( read2 ):
                    valid = False
                if fend2 != -self.find_fend( read1 ):
                    valid = False
            else:
                if fend1 != -self.find_fend( read2 ):
                    valid = False
        else:
            fend2, read2 = self.find_paired_fends( read2 )
            if fend2 != -self.find_fend( read1 ):
                valid = False
        if valid:
            self.valid_chimeric += 1
            self.valid += 1
            if read1 < read2:
                print >> self.output, '\t'.join(read1[:3] + read2[:3])
            else:
                print >> self.output, '\t'.join(read2[:3] + read1[:3])
        else:
            self.invalid_chimeric += 1

    def find_paired_fends( self, read ):
        cigar1 = self.parse_full_cigar( read[3] )
        cigar2 = self.parse_full_cigar( read[7] )
        if cigar1[0] < cigar2[0] and cigar1[1] > cigar2[1]:
            primary = read[:3]
            secondary = [read[4].strip('chr'), int(read[5]), read[6]]
        else:
            primary = read[4:7]
            secondary = [read[0].strip('chr'), int(read[1]), read[2]]
        if primary[2] == '-':
            primary[1] = str(int(primary[1]) + self.seq_len)
        if secondary[2] == '+':
            secondary[1] += self.seq_len
        fend = numpy.searchsorted(self.fends[secondary[0]], secondary[1])
        if secondary[2] == '-':
            fend = -fend
        return fend, primary

    def parse_full_cigar( self, cigar ):
        cigar = cigar.split('M')
        cigar[0] = cigar[0].rstrip('0123456789')
        left = 0
        index = cigar[0].find('H')
        if index == -1:
            index = cigar[0].find('S')
        while index != -1:
            left += int(cigar[0][:index])
            cigar[0] = cigar[0][(index + 1):]
            index = cigar[0].find('H')
            if index == -1:
                index = cigar[0].find('S')
        right = 0
        index = cigar[-1].find('H')
        if index == -1:
            index = cigar[-1].find('S')
        while index != -1:
            right += int(cigar[-1][:index])
            cigar[-1] = cigar[-1][(index + 1):]
            index = cigar[-1].find('H')
            if index == -1:
                index = cigar[-1].find('S')
        return [left, right]

    def find_fend( self, read ):
        if read[2] == '+':
            return numpy.searchsorted(self.fends[read[0].strip('chr')], int(read[1]))
        else:
            return -numpy.searchsorted(self.fends[read[0].strip('chr')], int(read[1]))

class Reader( object ):

    def __init__( self, fname ):
        self.fname = os.path.abspath( fname )
        self.reader = open( self.fname, 'r' )
        self.reads = {}
        self.strand = { '0':'+', '16':'-', '2048':'+', '2064':'-' }
        self.get_reads()
        self.reader.close()

    def get_reads( self ):
        for line in self.reader:
            split = line.rstrip('\n').split('\t')
            if split[0] in self.reads:
                temp = self.reads[split[0]]
                if not temp is None:
                    if len(temp) > 4:
                        self.reads[split[0]] = None
                    else:
                        self.reads[split[0]] = temp + [split[2], split[3], self.strand[split[1]], split[5]]
            else:
                self.reads[split[0]] = [split[2], split[3], self.strand[split[1]], split[5]]

if __name__ == "__main__":
    main()
