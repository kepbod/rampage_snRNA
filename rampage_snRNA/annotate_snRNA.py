#!/usr/bin/env python

'''
Usage: annotate_snRNA.py [options] -a snRNA <rampagedir>

Options:
    -h --help                      Show help message.
    --version                      Show version.
    -a snRNA                       snRNA annotations (BED format).
    --extend length                snRNA extended length. [default: 50]
    --span span                    span cutoff. [default: 1000]
    --coverage coverage            Coverage cutoff. [default: 0.5]
    -o out                         Output file. [default: snRNA_peak.txt]
'''

import os.path
from pybedtools import BedTool
import tempfile
from collections import defaultdict

__author__ = 'Xiao-Ou Zhang <xiaoou.zhang@umassmed.edu>'
__version__ = '0.0.1'


def annotate(options):
    '''
    Annotate expressed snRNA elements
    '''
    # check options
    extend = int(options['--extend'])
    span = int(options['--span'])
    coverage_percentage = float(options['--coverage'])
    output = options['-o']
    # prepare tmp files
    temp_dir = tempfile.mkdtemp()
    snRNA = temp_dir + '/snRNA.bed'
    # parse snRNA annotations
    with open(snRNA, 'w') as out:
        for chrom, start, end, name, strand in parse_snRNA(options):
            if strand == '+':
                start -= extend
                if start < 0:
                    continue
            else:
                end += extend
            out.write('%s\t%d\t%d\t%s\t0\t%s\n' % (chrom, start, end, name,
                                                   strand))
    # fetch snRNA peaks
    rampage_bed = BedTool(os.path.join(options['<rampagedir>'],
                                       'rampage_peak_entropy.txt'))
    snRNA_bed = BedTool(snRNA)
    snRNA_peak = rampage_bed.intersect(snRNA_bed, s=True, wa=True, wb=True)
    # filter snRNA peaks
    peak_lst = defaultdict(list)
    for p in snRNA_peak:
        if not check_span(p, span):
            continue
        summit = int(p[6])
        snRNA_chr = p[13]
        snRNA_start = int(p[14])
        snRNA_end = int(p[15])
        snRNA_name = p[16]
        snRNA_strand = p[18]
        if snRNA_strand == '+':
            snRNA_start += extend
            coverage = snRNA_end - summit
        else:
            snRNA_end -= extend
            coverage = summit - snRNA_start
        snRNA_len = snRNA_end - snRNA_start
        if coverage < snRNA_len * coverage_percentage:
            continue
        snRNA_info = '%s\t%d\t%d\t%s\t0\t%s' % (snRNA_chr, snRNA_start,
                                                snRNA_end, snRNA_name,
                                                snRNA_strand)
        peak_info = '\t'.join(p[:13])
        peak_lst[peak_info].append(snRNA_info)
    snRNA_exp = {}
    snRNA_info = {}
    for peak_info in peak_lst:
        expression = float(peak_info.split()[9])
        if len(peak_lst[peak_info]) == 1:
            snRNA = peak_lst[peak_info][0]
            if snRNA in snRNA_exp:
                if expression > snRNA_exp[snRNA]:
                    snRNA_exp[snRNA] = expression
                    snRNA_info[snRNA] = peak_info
            else:
                snRNA_exp[snRNA] = expression
                snRNA_info[snRNA] = peak_info
        else:
            strand, summit = peak_info.split()[5:7]
            summit = int(summit)
            for n, snRNA in enumerate(peak_lst[peak_info]):
                start, end = snRNA.split()[1:3]
                snRNA_site = int(start) if strand == '+' else int(end)
                snRNA_d = abs(snRNA_site - summit)
                if n == 0:
                    final_snRNA = snRNA
                    final_d = snRNA_d
                else:
                    if snRNA_d < final_d:
                        final_snRNA = snRNA
                        final_d = snRNA_d
            if final_snRNA in snRNA_exp:
                if expression > snRNA_exp[final_snRNA]:
                    snRNA_exp[final_snRNA] = expression
                    snRNA_info[final_snRNA] = peak_info
            else:
                snRNA_exp[final_snRNA] = expression
                snRNA_info[final_snRNA] = peak_info
    with open(output, 'w') as out:
        for snRNA in snRNA_info:
            out.write(snRNA_info[snRNA] + '\t' + snRNA + '\n')


def parse_snRNA(options):
    snRNA_f = options['-a']
    index = {'c': 0, 's': 1, 'e': 2, 'n': 3, 't': 5}
    with open(snRNA_f, 'r') as f:
        for line in f:
            items = line.rstrip().split()
            chrom = items[index['c']]
            start = int(items[index['s']])
            end = int(items[index['e']])
            name = items[index['n']]
            strand = items[index['t']]
            yield (chrom, start, end, name, strand)


def check_span(peak, span_cutoff):
    span = [int(x) for x in peak[11].split('|')]
    count = [float(x) for x in peak[12].split('|')]
    span_lst = {span[i]: count[i] for i in range(len(span))}
    effective_span = fetch_span(span_lst)
    if effective_span <= span_cutoff:
        return True
    else:
        return False


def fetch_span(span_lst):
    total = 0.75 * sum(span_lst.values())
    counts = 0
    span = None
    for s in sorted(span_lst):
        counts += span_lst[s]
        if counts >= total:
            span = s
            break
    return span


if __name__ == '__main__':
    from docopt import docopt
    annotate(docopt(__doc__, version=__version__))
