#!/usr/bin/env python
# coding=utf-8

import os
import sys
import time
import re
from itertools import chain
import argparse
import gzip
import traceback
import bz2file


def const():
    class _const:
        class ConstError(TypeError):
            pass

        class ConstCaseError(ConstError):
            pass

        def __setattr__(self, name, value):
            if self.__dict__.has_key(name):
                raise Exception(self.ConstError, "Can't change const.%s" % name)
            if not name.isupper():
                raise Exception(self.ConstCaseError, 'const name "%s" is not all uppercase!' % name)
            self.__dict__[name] = value

    sys.modules[__name__] = _const()


def readfasta(parser, path):
    out_data = []
    contig_seq = []
    contig_name = ''
    if path.endswith(".gz"):
        with gzip.open(path) as f:
            try:
                f.read(10)
            except IOError:
                print("Error in gzip file!")
                sys.exit(1)
            else:
                f.seek(0)
            for line in f:
                if line.startswith('>'):
                    contig_seqs = ''.join(contig_seq)
                    out_data.append([contig_name, contig_seqs])
                    contig_seq = []
                    contig_name = line.strip()
                else:
                    if parser.unzip:
                        contig_seq.append(line.strip())
                    else:
                        contig_seq.append(line.strip().upper())
    elif path.endswith(".bz2"):
        with bz2file.open(path) as f:
            try:
                f.read(10)
            except IOError:
                print("Error in bzip2 file!")
                sys.exit(1)
            else:
                f.seek(0)
            for line in f:
                if line.startswith('>'):
                    contig_seqs = ''.join(contig_seq)
                    out_data.append([contig_name, contig_seqs])
                    contig_seq = []
                    contig_name = line.strip()
                else:
                    if parser.unzip:
                        contig_seq.append(line.strip())
                    else:
                        contig_seq.append(line.strip().upper())
    else:
        with open(path) as f:
            for line in f:
                if line.startswith('>'):
                    contig_seqs = ''.join(contig_seq)
                    out_data.append([contig_name, contig_seqs])
                    contig_seq = []
                    contig_name = line.strip()
                else:
                    if parser.unzip:
                        contig_seq.append(line.strip())
                    else:
                        contig_seq.append(line.strip().upper())

    contig_seqs = ''.join(contig_seq)
    out_data.append([contig_name, contig_seqs])
    out_data.pop(0)
    return out_data

def prepare_argparesr():
    description = "%(prog)s -- Compressing sequencing data"
    usage = """ usage: %(prog)s <-f file> <-d> [-o output prefix]
    Example: %(prog)s -f hg19.fasta -o hg19
    """
    parser = argparse.ArgumentParser('DNACompressor', description=description, usage=usage)
    group = parser.add_mutually_exclusive_group()
    parser.add_argument('-f', '--file', help='input filename', action='store',
                        dest='input_file', required=True)
    parser.add_argument('-o', '--output', help='Output filename prefix(Default:noname)', action='store',
                        dest='output_file', default='noname')
    parser.add_argument('-d', '--decode', help='extract from the compressed file',
                        action='store_true', dest='unzip', default=False)
    parser.add_argument('-l', '--level', help='Compression level. Default:6', action='store',
                        default=6, type=int)
    return parser


def arg_validate(parser):
    parser = parser.parse_args()
    return parser


def translate():
    # first we have our score table
    list_3 = []
    list_2 = []
    list_1 = []
    comp_table = {}
    decomp_table = {}
    bases = ['T', 'C', 'A', 'G']
    for i in bases:
        list_1.append(i)
        for j in bases:
            list_2.append(''.join([i, j]))
            for k in bases:
                list_3.append(''.join([i, j, k]))

    for i1, i2 in enumerate(list_1):
        comp_table[i2] = chr(32 + i1)

    for j1, j2 in enumerate(list_2):
        comp_table[j2] = chr(48 + j1)

    for k1, k2 in enumerate(list_3):
        comp_table[k2] = chr(64 + k1)
    comp_table[''] = ''

    for key, value in comp_table.items():
        decomp_table[value] = key
    return comp_table, decomp_table


def make_Ntab(sequence):
    # findall 'N' and split other nucleotides

    N_seqs = re.findall(r'N+', sequence)
    N_seq = []
    for i in N_seqs:
        marks = '[N,' + str(len(i)) + ']'
        N_seq.append(marks)

    # Add one item for zip N_seq and A_seq
    N_seq.append('')
    A_seq = re.split(r'N+', sequence)

    if len(A_seq) == len(N_seq):
        pass
    else:
        print("Error in Tab length!")
        exit(0)
    return A_seq, N_seq


def find_zipN(sequence):
    zipN_seqs = re.findall(r'\[N\,\d{1,10}\]', sequence)
    N = 'N'
    zipN_seq = []
    for i in zipN_seqs:
        keys = eval(i)
        mulit_N = str(keys[0] * keys[1])
        zipN_seq.append(mulit_N)

    # for chain zipN_seq and zipA_Seq
    zipN_seq.append('')
    zipA_seq = re.split(r'\[N\,\d{1,10}\]', sequence)
    return zipA_seq, zipN_seq


def compress(files, comp_table):
    for sections in files:
        primary_sequence = sections[1]
        A_seq, N_seq = make_Ntab(primary_sequence)

        tab_seq = []
        for tabs in A_seq:
            seq = [comp_table[tabs[i:i + 3]] for i in range(0, len(tabs), 3)]
            tab_seq.append(''.join(seq))

        sections[1] = ''.join(list(chain.from_iterable(zip(tab_seq, N_seq))))
    return files


def decompress(files, decomp_table):
    for sections in files:
        zip_sequence = sections[1]
        zipA_seq, zipN_seq = find_zipN(zip_sequence)

        seq = []
        for keys in zipA_seq:
            seqs = [decomp_table[i] for i in keys]
            seq.append(''.join(seqs))

        sections[1] = ''.join(list(chain.from_iterable(zip(seq, zipN_seq))))
    return files


def run():
    parser = arg_validate(prepare_argparesr())
    inputfile = readfasta(parser, parser.input_file)

    comp_table, decomp_table = translate()

    if not parser.unzip:
        comp_data = compress(inputfile, comp_table)
    else:
        comp_data = decompress(inputfile, decomp_table)

    if not parser.unzip:
        f = open('%s.bc' % parser.output_file, 'w')
    else:
        f = open('%s.fasta' % parser.output_file, 'w')

    for items in comp_data:
        print >> f, items[0]
        temp = items[1]
        if parser.unzip:
            list_temp = [temp[i:i + 60] for i in range(0, len(temp), 60)]
            for key in list_temp:
                print >> f, key
        else:
            print >> f, temp

    f.close()

    if not parser.unzip:
        if os.path.exists("{0}.bc.bz2".format(parser.output_file)):
            sh("rm {0}.bc.bz2".format(parser.output_file))
        sh("bzip2 -{1} {0}.bc".format(parser.output_file, parser.level))

if __name__ == '__main__':
    sh = os.system
    begin_time = time.clock()
    start = time.strftime("%Y-%m-%d %X", time.localtime())

    try:
        run()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) Bye!\n")
        traceback.print_exc()
        sys.exit(0)

    end_time = time.clock()
    lapsed_time = end_time - begin_time

    print("\033[1;31;38m")
    print("#" * 50)
    print('Starts at \t' + start + '...')
    print('Ends at \t' + time.strftime("%Y-%m-%d %X", time.localtime()) + '...')
    print("Totally \t%.03f seconds lapsed" % lapsed_time)
    print("#" * 50)
    print("\033[0m")
