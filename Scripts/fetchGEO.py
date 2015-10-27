#!/usr/bin/env python

import sys
import os

import urllib2
from lxml import html

def main():
    gse_fname, outdir = sys.argv[1:3]
    REs = ['HindIII', 'NcoI', 'DpnII', 'MboI', 'MluI', 'BglII']
    datasets = {}
    for line in open(gse_fname):
        acc = line.strip('\n')
        if acc[0] == '#':
            continue
        if os.path.exists( "%s/%s.txt" % (outdir, acc) ):
            continue
        print acc
        datasets[acc] = load_samples(acc, REs, outdir)

def load_samples(acc, REs, outdir):
    results = {}
    samples = {}
    response = urllib2.urlopen("http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s" % acc)
    tree = html.fromstring(response.read())
    for e in tree.iter():
        if not e.text is None and e.text.count("Samples") > 0 and e.tag == 'td':
            c = e.getparent().getchildren()[1].getchildren()[0]
            for f in c.getchildren():
                g = f.getchildren()
                sra = g[0].getchildren()[0].text
                name = g[1].text
                name = name.replace(' ', '_').replace('(','_').replace(')','_').replace('__','_').replace('rep','Rep').replace('Replicate','Rep').replace(',','_')
                samples[sra] = { 'name': name }
            if len(e.getparent().getchildren()[1].getchildren()) >= 3:
                c = e.getparent().getchildren()[1].getchildren()[2].getchildren()[0]
                for f in c.getchildren():
                    g = f.getchildren()
                    sra = g[0].getchildren()[0].text
                    name = g[1].text
                    name = name.replace(' ', '_').replace('(','_').replace(')','_').replace('__','_').replace('rep','Rep').replace('Replicate','Rep').replace(',','_')
                    samples[sra] = { 'name': name } 
        elif not e.text is None and e.text.count("Contributor") > 0:
            c = e.getparent().getchildren()[1]
            authors = []
            for f in c.getchildren():
                authors.append(f.text)
            results['contributors'] = authors
        elif not e.text is None and e.text.count("Citation") > 0:
            c = e.getparent().getchildren()[1].getchildren()[0]
            if 'id' in c.keys():
                results['PMID'] = c.attrib['id']
    for sra in samples:
        response = urllib2.urlopen("http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s" % sra)
        tree = html.fromstring(response.read())
        for e in tree.iter():
            if not e.text is None and e.text == "SRA" and e.tag == 'td' and e.getparent().getchildren()[0].text == 'SRA':
                c = e.getparent().getchildren()[1].getchildren()[0]
                samples[sra]['srx'] = c.text
            if not e.text is None and e.text == "Organism" and e.tag == 'td':
                c = e.getparent().getchildren()[1].getchildren()[0]
                samples[sra]['organism'] = c.text
            if not e.text is None and e.text == "Source name" and e.tag == 'td':
                c = e.getparent().getchildren()[1]
                samples[sra]['source_name'] = c.text
            if not e.text is None and e.text == "Extraction protocol" and e.tag == 'td':
                c = e.getparent().getchildren()[1]
                res = []
                for r in REs:
                    if r in c.text:
                        res.append(r)
                for f in c.getchildren():
                    if f.text is None:
                        continue
                    for r in REs:
                        if r in c.text:
                            res.append(r)
                samples[sra]['REs'] = list(set(res))
            if not e.text is None and e.text == "Data processing" and e.tag == 'td':
                c = e.getparent().getchildren()[1]
                if 'REs' in samples[sra]:
                    res = samples[sra]['REs']
                else:
                    res = []
                for r in REs:
                    if r in c.text:
                        res.append(r)
                for f in c.getchildren():
                    if f.text is None:
                        continue
                    for r in REs:
                        if r in c.text:
                            res.append(r)
                samples[sra]['REs'] = list(set(res))
    for sra in samples:
        response = urllib2.urlopen("http://www.ncbi.nlm.nih.gov/sra?term=%s" % samples[sra]['srx'])
        tree = html.fromstring(response.read())
        srr = []
        for e in tree.iter():
            if 'class' in e.keys() and e.attrib['class'] == 'sra-run-list-header':
                for c in e.getparent().getparent().getchildren()[1].getchildren():
                    srr.append(c.getchildren()[0].getchildren()[0].text)
                samples[sra]['srr'] = srr
    results['samples'] = samples
    print_samples(acc, results, outdir)
    return results

def print_samples(acc, results, outdir):
    output = open("%s/%s.txt" % (outdir, acc), 'w')
    print >> output, "# Accession: %s" % acc
    contributors = ','.join(results['contributors'])
    if isinstance(contributors, str):
        contributors = unicode(contributors, 'utf-8')
    print >> output, "# Contributors: %s" % contributors.encode('utf-8')
    if 'PMID' in results:
        print >> output, "# Pubmed ID: %s" % results['PMID']
    samples = results['samples'].keys()
    samples.sort()
    for sample in samples:
        info = results['samples'][sample]
        if isinstance(info['name'], str):
            info['name'] = unicode(info['name'], 'utf-8')
        print >> output, "# Name: %s     Organism: %s     Source: %s     RE: %s" % (
            info['name'].encode('utf-8'), info['organism'], info['source_name'], ','.join(info['REs']))
        print >> output, "%s %s %s %s" % (sample, info['srx'], ','.join(info['srr']), ','.join(info['REs']))
    output.close()


if __name__ == "__main__":
    main()