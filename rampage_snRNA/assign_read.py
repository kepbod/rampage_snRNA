#!/usr/bin/env python

'''
Usage: assign_read.py [options] <rampagedir>

Options:
    -h --help                      Show help message.
    -v --version                   Show version.
    --tol TOL                      Termination for EM [default: 1e-9].
    --maxiter ITER                 Maximum iteration steps [default: 2000].
'''

import os.path
from pybedtools import BedTool
from seqlib.path import check_dir
from seqlib.em import squarem
from collections import defaultdict
import numpy as np

__author__ = 'Xiao-Ou Zhang <xiaoou.zhang@umassmed.edu>'
__version__ = '1.0.0'


def assign_read(options):
    '''
    Assign multiple mapping reads using EM algorithm
    '''
    # parse options
    folder = check_dir(options['<rampagedir>'])
    tol = float(options['--tol'])
    maxiter = int(options['--maxiter'])
    link_f = os.path.join(folder, 'rampage_link_count.bed')
    entropy_f = os.path.join(folder, 'rampage_peak_entropy.txt')
    with open(link_f, 'w') as link_out, open(entropy_f, 'w') as entropy_out:
        # overlap peaks with read paris
        link_bed = BedTool(os.path.join(folder, 'rampage_link.bed'))
        peak_bed = BedTool(os.path.join(folder, 'rampage_peaks.txt'))
        overlapped_peak = peak_bed.intersect(link_bed, s=True, wa=True,
                                             wb=True)
        # * peak_uniq_count: peak_info -> read_count
        peak_uniq_count = defaultdict(int)
        # * peak_multi_count: read_name -> {peak_info...}
        peak_multi_count = defaultdict(set)
        multi_read_info = {}  # for recording all the multiple read pairs
        peak_set = set()  # for recoding all the peak
        peak_span = {}
        for p in overlapped_peak:
            # peak info
            peak_chrom = p[0]
            peak_start = int(p[1])
            peak_end = int(p[2])
            peak_strand = p[5]
            peak_summit = int(p[6])
            peak_new_start = p[10]
            peak_new_end = p[11]
            peak_info = '\t'.join([peak_chrom, peak_new_start, peak_new_end,
                                   'peak\t0', peak_strand, str(peak_summit),
                                   str(peak_start), str(peak_end)])
            if peak_info not in peak_set:
                peak_set.add(peak_info)
                peak_span[peak_info] = defaultdict(float)
            # read pair info
            read_start = int(p[13])
            read_end = int(p[14])
            # check if read pair in peak region
            if peak_strand == '+':
                if read_start < peak_start or read_start > peak_end:
                    continue
            else:
                if read_end < peak_start or read_end > peak_end:
                    continue
            read_info = '\t'.join(p[12:])
            read_name = p[15]
            read_count = int(p[16])
            if peak_strand == '+':
                read_span = str(read_end - peak_summit)
            else:
                read_span = str(peak_summit - read_start)
            if read_name == 'unique':  # unique
                # write unique pairs
                link_out.write('{}\tunique\t{}\n'.format(read_info,
                                                         read_count))
                peak_uniq_count[peak_info] += read_count
                peak_span[peak_info][read_span] += read_count
            else:  # multiple
                peak_multi_count[read_name].add(peak_info)
                read_peak_id = '\t'.join([read_name, peak_info])
                # ! Note: read_span added here
                read_info += '\t{}'.format(read_span)
                multi_read_info[read_peak_id] = read_info
        # remove unique reads in multiple set
        for read_peak in multi_read_info:
            read, peak = read_peak_id.split('\t', 1)
            if len(peak_multi_count[read]) == 1:  # overlap only one peak
                read_info, read_span = multi_read_info[read_peak].rsplit('\t',
                                                                         1)
                link_out.write('{}\tunique\t1\n'.format(read_info))
                peak_span[peak][read_span] += 1
                del peak_multi_count[read]
        # * peak_lst: peak_info -> index
        peak_lst = {x: n for n, x in enumerate(peak_set)}  # make peak list
        # initiate expression
        peak_exp_init = update_exp(None, peak_lst, peak_uniq_count,
                                   peak_multi_count)
        # EM
        peak_exp, iter_num = squarem(peak_exp_init, update_exp, peak_lst,
                                     peak_uniq_count, peak_multi_count,
                                     tol=tol, maxiter=maxiter, verbose=True)
        print('EM converge with {} iterations'.format(iter_num))
        for read_peak, peak, exp in assign_multi_read(peak_exp, peak_lst,
                                                      peak_multi_count):
            read_info, read_span = multi_read_info[read_peak].rsplit('\t', 1)
            link_out.write('{}\tmultiple\t{}\n'.format(read_info, exp))
            peak_span[peak][read_span] += exp
        # calculate entropy
        for peak in peak_span:
            if peak_span[peak]:
                entropy = 0
                total = sum(peak_span[peak].values())
                for span in peak_span[peak]:
                    p = peak_span[peak][span] / total
                    entropy += p * np.log2(p)
                entropy = -entropy
                out = '\t'.join([peak, str(total), str(entropy),
                                 '|'.join(peak_span[peak]),
                                 '|'.join(str(x)
                                          for x in peak_span[peak].values())])
                entropy_out.write('{}\n'.format(out))


def update_exp(peak_exp_pre, peak_lst, uniq, multi):
    if peak_exp_pre is None:
        init_flag = True
    else:
        init_flag = False
    peak_exp = np.zeros(len(peak_lst))
    peak_multi_exp = defaultdict(float)
    # add unique reads first
    for peak in uniq:
        peak_exp[peak_lst[peak]] += uniq[peak]
    # then distribute multiple reads
    for read in multi:
        if init_flag:
            fetch_exp_init(peak_multi_exp, multi[read], peak_exp, peak_lst)
        else:
            fetch_exp(peak_multi_exp, multi[read], peak_exp_pre, peak_lst)
    for peak in peak_multi_exp:
        peak_exp[peak_lst[peak]] += peak_multi_exp[peak]
    peak_exp /= peak_exp.sum()
    return peak_exp


def fetch_exp_init(peak_multi_exp, peak_set, peak_exp, peak_lst):
    number = len(peak_set)
    total = 1
    peak_with_uniq = {}
    for peak in peak_set:
        exp = peak_exp[peak_lst[peak]]
        if exp == 0:  # no unique reads
            peak_multi_exp[peak] += 1 / number
            total -= 1 / number
        else:  # have unique reads
            peak_with_uniq[peak] = exp
    peak_with_uniq_total = sum(peak_with_uniq.values())
    for peak in peak_with_uniq:
        peak_multi_exp[peak] += (total * peak_with_uniq[peak] /
                                 peak_with_uniq_total)


def fetch_exp(peak_multi_exp, peak_set, peak_exp, peak_lst):
    peak_set_exp = {peak: peak_exp[peak_lst[peak]] for peak in peak_set}
    peak_set_exp_total = sum(peak_set_exp.values())
    for peak in peak_set_exp:
        if peak_set_exp_total == 0:
            peak_multi_exp[peak] += 1 / len(peak_set)
        else:
            peak_multi_exp[peak] += peak_set_exp[peak] / peak_set_exp_total


def assign_multi_read(peak_exp, peak_lst, multi):
    for read in multi:
        peak_set_exp = {peak: peak_exp[peak_lst[peak]] for peak in multi[read]}
        peak_set_exp_total = sum(peak_set_exp.values())
        for peak in peak_set_exp:
            exp = peak_set_exp[peak] / peak_set_exp_total
            if exp == 0:
                continue
            read_peak = '\t'.join([read, peak])
            yield read_peak, peak, exp


if __name__ == "__main__":
    from docopt import docopt
    assign_read(docopt(__doc__, version=__version__))
