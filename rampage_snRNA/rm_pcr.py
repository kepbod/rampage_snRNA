#!/usr/bin/env python

'''
Usage: rm_pcr.py [options] <rampage>...

Options:
    -h --help                      Show help message.
    --version                      Show version.
    -a snRNA                       snRNA annotations (BED format).
    -p THREAD --thread=THREAD      Threads. [default: 5]
    -o OUTPUT --output=OUTPUT      Output directory. [default: rampage_peak]
'''

import os
import os.path
from collections import defaultdict
from seqlib.path import create_dir
from seqlib.ngs import check_bam
from seqlib.seq import dna_to_rna
from seqlib.interval import Interval

__author__ = 'Xiao-Ou Zhang <xiaoou.zhang@umassmed.edu>'
__version__ = '1.0.0'


def rm_pcr(options):
    '''
    Remove PCR duplicates
    '''
    # parse options
    rampage_lst = options['<rampage>']
    thread = int(options['--thread'])
    folder = create_dir(options['--output'])
    # fetch promoter regions of snRNAs
    snRNA_promoter = fetch_promoter(options['-a'])
    # remove PCR duplicates
    from multiprocessing import Pool
    p = Pool(thread)
    result = []
    for rampage in rampage_lst:  # for each rampage replicate
        for chrom in check_bam(rampage).references:  # for each chromosome
            if chrom.startswith('chr'):  # only parse normal chromosome
                promoter = {}
                promoter['+'] = snRNA_promoter['|'.join([chrom, '+'])]
                promoter['-'] = snRNA_promoter['|'.join([chrom, '-'])]
                result.append(p.apply_async(remove_pcr, args=(rampage,
                                                              promoter,
                                                              chrom,)))
    p.close()
    p.join()
    uniq_pairs = []
    multi_pairs = defaultdict(set)
    for r in result:
        uniq, multi = r.get()
        uniq_pairs.extend(uniq)
        for name in multi:
            multi_pairs[name].update(multi[name])
    write_signal(uniq_pairs, multi_pairs, folder)


def fetch_promoter(bed):
    promoter_region = defaultdict(list)
    with open(bed, 'r') as bed_f:
        for line in bed_f:
            chrom, start, end, _, _, strand = line.rstrip().split()
            start = int(start)
            end = int(end)
            if strand == '+':
                s = start - 100
                e = start + 100
            else:
                s = end - 100
                e = end + 100
            promoter_region['|'.join([chrom, strand])].append([s, e])
    for chrom in promoter_region:
        promoter_region[chrom] = Interval(promoter_region[chrom])
    return promoter_region


def write_signal(uniq_pairs, multi_pairs, folder):
    bed5p = open(os.path.join(folder, 'rampage_plus_5end.bed'), 'w')
    bed5m = open(os.path.join(folder, 'rampage_minus_5end.bed'), 'w')
    link = open(os.path.join(folder, 'rampage_link.bed'), 'w')
    read_pairs = defaultdict(int)
    for pair in uniq_pairs:
        info = pair.rsplit('\t', 1)[0]  # remove barcode info
        info += '\tunique'
        read_pairs[info] += 1
    pos_list = []
    # remove PCR duplicates for multiple mapping reads
    for name in multi_pairs:
        pos = multi_pairs[name]
        if pos not in pos_list:  # not PCR duplicates
            pos_list.append(pos)
            for p in pos:
                info = p.rsplit('\t', 1)[0]  # remove barcode info
                info += '\t{}'.format(name)
                read_pairs[info] += 1
    for pair in read_pairs:
        read_count = read_pairs[pair]
        pair_info = pair.split()
        r1_chrom, r1_start, r1_end, strand = pair_info[:4]  # read1
        r2_chrom, r2_start, r2_end, name = pair_info[4:]  # read2
        if strand == '+':
            start = int(r1_start)
            end = int(r2_end)
            if end - start > 2000:  # exclude very long pairs
                continue
            bed5p.write('%s\t%d\t%d\t5end\t0\t%s\n' % (r1_chrom, start,
                                                       start + 1, strand) *
                        read_count)
            offset = '0,' + str(end - start - 1)
            link.write('\t'.join([r1_chrom, r1_start, r2_end, name,
                                 str(read_count), strand, r1_start, r1_start,
                                 '0,0,0', '2', '1,1', offset]) + '\n')
        else:
            start = int(r2_start)
            end = int(r1_end)
            if end - start > 2000:  # exclude very long pairs
                continue
            bed5m.write('%s\t%d\t%d\t5end\t0\t%s\n' % (r1_chrom, end - 1,
                                                       end, strand) *
                        read_count)
            offset = '0,' + str(end - start - 1)
            link.write('\t'.join([r1_chrom, r2_start, r1_end, name,
                                 str(read_count), strand, r2_start, r2_start,
                                 '0,0,0', '2', '1,1', offset]) + '\n')


def remove_pcr(bam_f, promoter, chrom):
    if type(bam_f) is list:  # multiple bam files
        uniq_pairs = []
        multi_pairs = defaultdict(set)
        for f in bam_f:
            bam = check_bam(f)
            read1 = fetch_read1(bam, chrom, promoter)  # fetch read1
            # fetch read2
            uniq, multi = fetch_read2(bam, read1, chrom)
            uniq_pairs.extend(uniq)
            for name in multi:
                multi_pairs[name].update(multi[name])
    else:  # single bam file
        bam = check_bam(bam_f)
        read1 = fetch_read1(bam, chrom, promoter)  # fetch read1
        uniq_pairs, multi_pairs = fetch_read2(bam, read1, chrom)  # fetch read2
    return uniq_pairs, multi_pairs


def fetch_read1(bam, chrom, promoter):
    read1_lst = {}
    for read in bam.fetch(chrom):
        # ! Note: we did not exclude secondary alignments
        # ! because we want to determine potential locations
        # ! for multiple mapping
        if read.is_read2:  # not read2
            continue
        if not read.is_proper_pair:  # not proper pair
            continue
        if read.get_tag('NH') == 1:  # unique mapping
            uniq_flag = True
        else:  # multiple mapping
            uniq_flag = False
        chrom = read.reference_name
        mate_chrom = read.next_reference_name
        if chrom != mate_chrom:  # not same chromosome
            continue
        if not read.is_reverse and read.mate_is_reverse:
            strand = '+'
            mate_strand = '-'
        elif read.is_reverse and not read.mate_is_reverse:
            strand = '-'
            mate_strand = '+'
        else:
            continue
        start = read.reference_start
        end = read.reference_end
        # check if in promoter regions
        if strand == '+':
            if [start, start + 1] not in promoter[strand]:
                continue
        else:
            if [end - 1, end] not in promoter[strand]:
                continue
        mate_pos = str(read.next_reference_start)
        name = read.query_name
        # * read_id: read_name, chrom, read2_pos, read2_strand, read1_pos
        read_id = '\t'.join([name, mate_chrom, mate_pos, mate_strand,
                             str(start)])
        read1_lst[read_id] = [chrom, str(start), str(end), strand, uniq_flag]
    return read1_lst


def fetch_read2(bam, read1, chrom):
    collapsed_pairs_uniq = set()
    collapsed_pairs_multi = defaultdict(set)
    for read in bam.fetch(chrom):
        # not read1
        if read.is_read1:
            continue
        if not read.is_proper_pair:  # not proper pair
            continue
        name = read.query_name
        chrom = read.reference_name
        start = str(read.reference_start)
        # parse strand info
        strand = '+' if not read.is_reverse else '-'
        mate_pos = str(read.next_reference_start)
        # * read_id: read_name, chrom, read2_pos, read2_strand, read1_pos
        read_id = '\t'.join([name, chrom, start, strand, mate_pos])
        if read_id not in read1:
            continue
        end = str(read.reference_end)
        if strand == '+':
            barcode = dna_to_rna(read.query_sequence[:15])
        else:
            barcode = dna_to_rna(read.query_sequence[-15:],
                                 strand=strand)
        # remove PCR duplicates
        # * pari_info: chrom, read1_start, read1_end, read1_strand
        # *            chrom, read2_start, read2_end, barcode
        pair_info = '\t'.join(read1[read_id][:4] +
                              [chrom, start, end, barcode])
        if read1[read_id][4]:  # unique mapping
            collapsed_pairs_uniq.add(pair_info)
        else:  # multiple mapping
            collapsed_pairs_multi[name].add(pair_info)
    return collapsed_pairs_uniq, collapsed_pairs_multi


if __name__ == '__main__':
    from docopt import docopt
    rm_pcr(docopt(__doc__, version=__version__))
