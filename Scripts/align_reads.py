#!/usr/bin/env python

import sys
import os
import subprocess
import argparse

def main():
    parser = generate_parser()
    args = parser.parse_args()
    aligner = Aligner( args.input, args.output, args.index, args.threads, args.buffer )
    return None

def generate_parser():
    """Generate an argument parser."""
    description = "%(prog)s -- Create a raw file of paired aligned reads for a HiC experiment from fastq or sam files."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument( "-i", "--index", dest="index", required=False, type=str, default=None,
        action='store', help="BWA index if using unaligned fastq files. [default: %(default)s]" )
    parser.add_argument( "-t", "--threads", dest="threads", required=False, type=int, default=2,
        action='store', help="The number of total threads to use if aligning reads. [default: %(default)s]" )
    parser.add_argument( "-b", "--buffer", dest="buffer", required=False, type=int, default=4000000,
        action='store', help="The number of lines to be placed in a temporary file for mapping. [default: %(default)s]" )
    parser.add_argument( dest="input", type=str, action='store', help="Fastq file." )
    parser.add_argument( dest="output", type=str, action='store', help="Sam file prefix." )
    return parser

class Aligner( object ):

    def __init__( self, fname, outfname, index=None, num_threads=1, buffersize=4000000 ):
        self.fname = os.path.abspath( fname )
        self.index = os.path.abspath( index )
        self.outfname = outfname
        self.bufferfname = '%s.temp' % self.fname 
        self.buffersize = (buffersize / 4) * 4
        self.buffernum = 0
        self.fastq = open( self.fname, 'r' )
        self.fastq_line = self.fastq.readline()
        self.seq_len = self.get_seq_len()
        self.num_threads = max( num_threads - 1, 1 )
        self.strand = { '0': '+', '16': '-', '2048': '+', '2064': '-' }
        self.primary = { '0': None, '16': None }
        self.secondary = { '2048': None, '2064': None }
        self.total = 0
        self.aligned = 0
        self.unfiltered = 0
        self.align_reads()

    def get_seq_len( self ):
        try:
            seq_len = int(self.fastq_line.rstrip('\n').split('length=')[1])
        except:
            temp_in = open(self.fname)
            temp_in.readline()
            seq_len = len(temp_in.readline().rstrip('\n'))
            temp_in.close()
        return seq_len

    def align_reads( self ):
        while self.fill_buffer():
            self.bwa = subprocess.Popen( ['bwa', 'mem', '-c 1', '-t',
                                         str(self.num_threads), self.index, self.bufferfname],
                                         stdout=subprocess.PIPE, stderr=subprocess.PIPE )
            self.reader = self.bwa.stdout
            self.read_header()
            self.filter_reads()
            self.output.close()
        self.write_stats()
        return 0

    def fill_buffer( self ):
        output = open( self.bufferfname, 'w' )
        i = 0
        while self.fastq_line and i < self.buffersize:
            output.write( self.fastq_line )
            self.fastq_line = self.fastq.readline()
            i += 1
        output.close()
        if i > 0:
            self.output = open( "%s.%i" % (self.outfname, self.buffernum), 'w' )
            self.buffernum += 1
            return True
        else:
            return False

    def read_header( self ):
        self.line = self.reader.readline()
        while self.line and self.line[0] == '@':
            self.line = self.reader.readline()
        return None

    def filter_reads( self ):
        while self.line:
            split = self.line.split('\t')
            alignment = split[1]
            if alignment not in self.secondary:
                self.total += 1
            if alignment in self.strand:
                if alignment in self.primary:
                    self.aligned += 1
                if split[4] != '0':
                    print >> self.output, '\t'.join(split[:6] + split[8:9])
                    self.unfiltered += 1
            self.line = self.reader.readline()
        return False

    def write_stats( self ):
        print >> sys.stderr, ("Sequence length: %i\n") % (self.seq_len),
        print >> sys.stderr, ("Total reads: %i\n") % (self.total),
        print >> sys.stderr, ("Aligned reads: %i (%0.2f%%)\n") % (self.aligned, self.aligned * 100.0 / float(self.total)),
        print >> sys.stderr, ("High-quality reads: %i (%0.2f%%)\n") % (self.unfiltered, self.unfiltered * 100.0 / float(self.total)),


if __name__ == "__main__":
    main()
