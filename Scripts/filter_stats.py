#!/usr/bin/env python

import sys
import os

import numpy

import hifive

def main():
    in_prefix, out_fname = sys.argv[1:3]
    data_fname = "%s.hcd" % in_prefix
    project_fname = "%s.hcp" % in_prefix
    if not os.path.exists(project_fname):
        hic = hifive.HiC(project_fname, 'w')
        hic.load_data(data_fname)
        hic.save()
    else:
        hic = hifive.HiC(project_fname)
    hic.filter.fill(1)
    results = [hic.filter.shape[0]]
    cont = True
    i = 1
    while cont:
        hic.filter_fends(mininteractions=i, mindistance=0, maxdistance=0)
        results.append(numpy.sum(hic.filter))
        i += 1
        if results[-1] == 0:
            cont = False
    numpy.savetxt(out_fname, numpy.array(results, dtype=numpy.int32))


if __name__ == "__main__":
    main()