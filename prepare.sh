#!/bin/bash

MAX_PLASMID_LEN=500000
RANDOM_SEED=42
CLASS_SAMPLES=0

if [ ! -e plasmid_seq.tsv ]; then
	./plasmid_seq_extractor.py --plasmids_fna plasmids/ncbi_plasmid.fna --plasmids_gbff plasmids/ncbi_plasmid.gbff \
			--output_file plasmid_seq.fna --max_length $MAX_PLASMID_LEN
fi
