#!/usr/bin/env python

import sys
import subprocess


def main():
    in_prefix, out_dir = sys.argv[1:3]
    stats = find_partitions(in_prefix)
    reads = {}
    num_reads = 0
    for i in range(stats['partitions']):
        fname = "%s.%i" % (in_prefix, i)
        for line in open(fname):
            num_reads += 1
            reads[tuple(line.rstrip('\n').split('\t'))] = None
    duplicates = num_reads - len(reads)
    for key in reads.keys():
        print '\t'.join(key)
    total = stats['total']
    aligned1 = stats['aligned1']
    aligned2 = stats['aligned2']
    quality1 = stats['quality1']
    quality2 = stats['quality2']
    chimeric = 0
    vchimeric = 0
    paired = 0
    valid = 0
    for i in range(stats['partitions']):
        fname = "%s.%i" % (in_prefix.replace('raw', 'stats'), i)
        for line in open(fname):
            temp = line.rstrip('\n').split(': ')
            if temp[0] == 'Paired reads':
                paired += int(temp[1])
            elif temp[0] == 'Valid paired reads':
                valid += int(temp[1])
            elif temp[0] == 'Total chimeric reads':
                chimeric += int(temp[1])
            elif temp[0] == 'Valid chimeric reads':
                vchimeric += int(temp[1])
    valid -= duplicates
    print >> sys.stderr, ("Total reads: %i\n") % (total),
    print >> sys.stderr, ("First pair aligned: %i (%0.2f%%)\n") % (aligned1, aligned1 * 100.0 / float(total)),
    print >> sys.stderr, ("Second pair aligned: %i (%0.2f%%)\n") % (aligned2, aligned2 * 100.0 / float(total)),
    print >> sys.stderr, ("First pair quality aligned: %i (%0.2f%%)\n") % (quality1, quality1 * 100.0 / float(total)),
    print >> sys.stderr, ("Second pair quality aligned: %i (%0.2f%%)\n") % (quality2, quality2 * 100.0 / float(total)),
    print >> sys.stderr, ("PCR duplicate paired reads: %i (%0.2f%%)\n") % (duplicates, duplicates * 100.0 / float(total)),
    print >> sys.stderr, ("Chimeric paired reads: %i (%0.2f%%)\n") % (chimeric, chimeric * 100.0 / float(total)),
    print >> sys.stderr, ("Valid  chimeric paired reads: %i (%0.2f%%)\n") % (vchimeric, vchimeric * 100.0 / float(total)),
    print >> sys.stderr, ("Total mapped/paired reads: %i (%0.2f%%)\n") % (paired, paired * 100.0 / float(total)),
    print >> sys.stderr, ("Valid mapped/paired reads: %i (%0.2f%%)\n") % (valid, valid * 100.0 / float(total)),
    if isinstance(valid, int) and valid > 0:
        subprocess.call("mv %s %s/" % (in_prefix, out_dir.rstrip('/')), shell=True)
        subprocess.call("mv %s %s/" % (in_prefix.replace('raw', 'stats'), out_dir.rstrip('/')), shell=True)
        if valid / float(total) >= 0.6:
            clean_files(in_prefix.split('.raw')[0])
    return None

def clean_files(srr):
    subprocess.call("rm -f %s.*" % srr, shell=True)
    subprocess.call("rm -f %s_*" % srr.replace('Mapping', 'Fastq'), shell=True)
    subprocess.call("rm -f %s.*" % srr.replace('Mapping', 'tmp'), shell=True)

def find_partitions(prefix):
    stats = {}
    for line in open("%s_1.stats" % prefix.split('.raw')[0]):
        temp = line.rstrip('\n').split(': ')
        if temp[0] == 'Total reads':
            stats['total'] = int(temp[1].strip(' '))
            stats['partitions'] = (stats['total'] - 1) / 2000000 + 1
        elif temp[0] == 'Aligned reads':
            stats['aligned1'] = int(temp[1].strip(' ').split(' ')[0])
        elif temp[0] == 'High-quality reads':
            stats['quality1'] = int(temp[1].strip(' ').split(' ')[0])
    for line in open("%s_2.stats" % prefix.split('.raw')[0]):
        temp = line.rstrip('\n').split(': ')
        if temp[0] == 'Aligned reads':
            stats['aligned2'] = int(temp[1].strip(' ').split(' ')[0])
        elif temp[0] == 'High-quality reads':
            stats['quality2'] = int(temp[1].strip(' ').split(' ')[0])
    return stats

if __name__ == "__main__":
    main()
