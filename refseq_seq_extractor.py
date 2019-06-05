#!/usr/bin/env python

from datetime import datetime
from functools import partial
import io
import os
import shutil

from Bio import SeqIO
import click
import multiprocessing
import numpy as np
import pandas as pd


def process_pool_init(lock, outfile, include_all):
    global output_file_lock
    output_file_lock = lock
    global output_f
    output_f = outfile
    global include_all_loci
    include_all_loci = include_all

def gff_seq_extract(filepath, min_length, max_length):
    output = io.StringIO()
    try:
        for seq_record in SeqIO.parse(filepath, 'genbank'):
            taxonomy = seq_record.annotations['taxonomy']
            if not include_all_loci and \
                    ('plasmid' in seq_record.description or
                        'extrachromosomal' in seq_record.description):
                continue
            if max_length == 0:
                less_than_max = 1
            else:
                less_than_max = 0
            if max_length > 0 and len(seq_record.seq) <= max_length:
                less_than_max = 1
            if min_length == 0:
                greater_than_min = 1
            else:
                greater_than_min = 0
            if min_length > 0 and len(seq_record.seq) >= min_length:
                greater_than_min = 1
            if less_than_max and greater_than_min:
                output.write(f">{seq_record.id}\t{taxonomy[0] if len(taxonomy) > 0 else np.nan}\t{taxonomy[1] if len(taxonomy) > 1 else np.nan}\t{taxonomy[2] if len(taxonomy) > 2 else np.nan}\t{taxonomy[3] if len(taxonomy) > 3 else np.nan}\t{taxonomy[4] if len(taxonomy) > 4 else np.nan}\t{taxonomy[5] if len(taxonomy) > 5 else np.nan}\n{seq_record.seq}\n")
    except AttributeError:
        print(f'parsing of file {filepath} failed')

    output.seek(0)
    output_file_lock.acquire()
    shutil.copyfileobj(output, output_f)
    output_f.flush()
    output_file_lock.release()


@click.command()
@click.option('-t', '--threads', 'threads', default=16, type=int,
        help='number of parallel Genbank parser threads')
@click.option('-r', '--refseq_path', 'refseq_path', type=str, required=True,
        help='path to the root of the refseq download')
@click.option('-o', '--output_file', 'output_file', type=str, required=True,
        help='name of the output file w/ taxonomy annotations and sequences')
@click.option('-g', '--genbank_postfix', 'genbank_postfix', type=str,
        default='.gbff', show_default=True,
        help='postfix for Genbank files')
@click.option('-i', '--include_all', 'include_all', type=bool,
        default=False, show_default=True,
        help='should plasmids and extrachromosomal elements be included')
@click.option('-l', '--max_length', 'max_length', type=int, required=False,
        default=0, help='maximum length of sequences to save to output')
@click.option('-i', '--min_length', 'min_length', type=int, required=False,
        default=0, help='minimum length of sequences to save to output')
def refseq_seq_extractor(threads, refseq_path, output_file, genbank_postfix, include_all, min_length, max_length):
    """Extract the taxonomy and CDS sequences from a collection of GBFF files"""

    start_time = datetime.now()

    print(f'finding Genbank files in {refseq_path}')
    gbff_files = []
    for root, dirs, files in os.walk(refseq_path):
        for file_ in files:
            filepath = os.path.join(root, file_)
            if filepath.endswith(genbank_postfix):
                gbff_files.append(filepath)

    print(f'scanning {len(gbff_files)} files to extract taxonomy and CDS sequences')
    print(f'using {threads} parallel parsers')

    f = open(output_file, 'w')
    #f.write('gff_file\tid\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tproduct_id\tsequence\n')
    f.flush()

    lock = multiprocessing.Lock()
    process_pool = multiprocessing.Pool(threads,
            initializer=process_pool_init, initargs=(lock, f, include_all, ))
    gff_seq_extract_f = partial(gff_seq_extract, min_length=min_length, max_length=max_length)
    process_pool.map(gff_seq_extract_f, gbff_files)
    process_pool.close()
    process_pool.join()

    f.close()

    stop_time = datetime.now()
    total_time = stop_time - start_time
    print(f'run time was: {total_time}')


if __name__ == '__main__':
    refseq_seq_extractor()
