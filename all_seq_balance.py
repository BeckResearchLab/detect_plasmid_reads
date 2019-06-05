#!/usr/bin/env python

import click
import numpy as np
import pandas as pd


@click.command()
@click.option('-o', '--output_file', 'output_file', type=str, required=True,
        help='name of the output file containing balanced samples')
@click.option('-i', '--input_file', 'input_file', type=str, required=True,
        help='the location of the taxonomy annotated sequences')
@click.option('-r', '--random_seed', 'random_seed', type=int, default=42,
        show_default=True,
        help='random seed to be used for shuffling and sampling data partitions')
@click.option('-s', '--positive_samples', 'positive_samples', type=int,
        default=0, show_default=True,
        help='limit the number of positive samples to this number (disabled = 0)')
def all_seq_balance(input_file, output_file, random_seed, positive_samples):
    """Balance samples of sequences for a binary classification on the is_plasmid field"""

    print(f'reading seq tsv data file {input_file}')
    df = pd.read_csv(input_file, sep='\t').filter(items=['sequence', 'is_plasmid'])
    print(f'there were {df.shape[0]} sequences in {input_file}')

    print(f'finding plasmid sequences')
    positives = df.loc[df['is_plasmid'] == 1]
    print(f'found {positives.shape[0]} samples')
    if positive_samples > 0:
        print(f'down sampling positive samples to {positive_samples} while shuffling')
        positives = positives.sample(n=positive_samples, random_state=random_seed)
    print(f'finding non-plasmid sequences')
    negatives = df.loc[df['is_plasmid'] == 0]
    print(f'found {negatives.shape[0]} samples to be randomly sampled down to {positives.shape[0]} samples')
    negatives = negatives.sample(n=positives.shape[0], random_state=random_seed)
    print(f'concatenating positive and negative samples')
    df = positives.append(negatives)
    # attempt to allow garbage collection before shuffling
    positives = None
    negatives = None
    print(f'shuffling the order of samples')
    df = df.sample(frac=1, random_state=random_seed).reset_index(drop=True)
    print(f'saving {df.shape[0]} balanced samples to {output_file}')
    df.to_csv(output_file, sep='\t', index=False)


if __name__ == '__main__':
    all_seq_balance()
