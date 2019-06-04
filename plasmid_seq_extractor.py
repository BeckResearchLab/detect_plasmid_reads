#!/usr/bin/env python

from datetime import datetime
import io
import os
import shutil

from Bio import SeqIO
import click
import numpy as np
import pandas as pd


@click.command()
@click.option('-f', '--plasmids_fna', 'plasmids_fna_path', type=str, required=True,
        help='path to the FASTA nucleotide file of the plasmids')
@click.option('-g', '--plasmids_gbff', 'plasmids_gbff_path', type=str,
        help='path to the Genbank file of the plasmids')
@click.option('-o', '--output_file', 'output_file', type=str, required=True,
        help='name of the output file w/ plasmid sequences')
@click.option('-l', '--max_length', 'max_length', type=int, required=False,
        default=0, help='maximum length of sequences to save to output')
@click.option('-i', '--min_length', 'min_length', type=int, required=False,
        default=0, help='minimum length of sequences to save to output')
def plasmid_seq_extractor(plasmids_fna_path, plasmids_gbff_path, output_file, max_length, min_length):
    """Extract the taxonomy and CDS sequences from a GBFF and FNA file for plasmids"""

    if max_length > 0 and min_length > 0:
        if min_length > max_length:
            raise ValueError('incompatible min and max length', f"max length ({max_length}) can't be less than min length ({min_length})")

    start_time = datetime.now()

    print(f'writing records to {output_file}')
    f = open(output_file, 'w')
    f.write('gff_file\tid\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tproduct_id\tsequence\tis_plasmid\n')
    f.flush()

    print(f'reading FASTA file {plasmids_fna_path}')
    fna_dict = SeqIO.to_dict(SeqIO.parse(plasmids_fna_path, 'fasta'))

    print(f'reading GBFF file {plasmids_gbff_path}')
    try:
        for seq_record in SeqIO.parse(plasmids_gbff_path, 'genbank'):
            taxonomy = seq_record.annotations['taxonomy']
            if 'contig' in seq_record.annotations:
                seqid = 'ref|' + seq_record.id + '|'
                if seqid in fna_dict:
                    seq_record.seq = fna_dict[seqid].seq
                else:
                    print(f'sequence {seqid} is missing from FASTA, skipping this record')
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
                f.write(f">{seq_record.id}\t{taxonomy[0] if len(taxonomy) > 0 else np.nan}\t{taxonomy[1] if len(taxonomy) > 1 else np.nan}\t{taxonomy[2] if len(taxonomy) > 2 else np.nan}\t{taxonomy[3] if len(taxonomy) > 3 else np.nan}\t{taxonomy[4] if len(taxonomy) > 4 else np.nan}\t{taxonomy[5] if len(taxonomy) > 5 else np.nan}\t{feature.qualifiers['protein_id'][0] if 'protein_id' in feature.qualifiers else np.nan}\n{seq_record.seq}\n")
    except AttributeError:
        print(f'parsing of file {filepath} failed')

    f.close()

    stop_time = datetime.now()
    total_time = stop_time - start_time
    print(f'run time was: {total_time}')


if __name__ == '__main__':
    plasmid_seq_extractor()
